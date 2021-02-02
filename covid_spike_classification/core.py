"""Core functions for covid spike classification."""

import glob
import os
import shutil
import subprocess
import sys
import zipfile

from .config import CSCConfig

from Bio.Seq import Seq

REGIONS = {
    # The important ones the system asks for
    "N439K": "NC_045512:22877-22879",
    "Y453F": "NC_045512:22919-22921",
    "E484K": "NC_045512:23012-23014",
    "N501Y": "NC_045512:23063-23065",
    "P681H": "NC_045512:23603-23605",
    # Bonus mutations to help call variants
    "L452R": "NC_045512:22916-22918",
    "S477N": "NC_045512:22991-22993",
    "A570D": "NC_045512:23270-23272",
    "D614G": "NC_045512:23402-23404",
    "A626S": "NC_045512:23438-23441",
    "H655Y": "NC_045512:23525-23527",
    "I692V": "NC_045512:23636-23638",
    "A701V": "NC_045512:23663-23665",
    "T716I": "NC_045512:23708-23710",
}


KNOWN_VARIANTS = {
    "B.1.1.7": {
        "N501Y",
        "A570D",
        "P681H",
        "T716I",
    },
    "B.1.351": {
        "E484K",
        "N501Y",
        "A701V",
    },
    "P1": {
        "E484K",
        "N501Y",
        "H655Y",
    },
}


class PileupFailedError(RuntimeError):
    pass


class BaseDeletedError(RuntimeError):
    pass


def basecall(tmpdir, config):
    if config.input_format != "ab1":
        return
    fastq_dir = os.path.join(tmpdir, "fastqs")
    os.makedirs(fastq_dir)

    ab1_dir = _extract_if_zip(tmpdir, config)

    for sanger_file in glob.glob(os.path.join(ab1_dir, "*.ab1")):
        base_name = os.path.basename(sanger_file)
        cmd = ["tracy", "basecall", "-f", "fastq", "-o", os.path.join(fastq_dir, f"{base_name}.fastq"), sanger_file]
        kwargs = {}
        if config.quiet:
            kwargs["stdout"] = subprocess.DEVNULL
            kwargs["stderr"] = subprocess.DEVNULL
        subprocess.check_call(cmd, **kwargs)

    shutil.copytree(fastq_dir, config.outdir, dirs_exist_ok=True)


def map_reads(tmpdir, config):

    if config.input_format == "ab1":
        # fastqs live in the tmpdir
        sequence_dir = os.path.join(tmpdir, "fastqs")
        file_ending = "fastq"
    else:
        sequence_dir = _extract_if_zip(tmpdir, config)
        file_ending = config.input_format


    bam_dir = os.path.join(tmpdir, "bams")
    os.makedirs(bam_dir)

    # ditch the .fasta file ending
    name, _ = os.path.splitext(config.reference)
    ref = f"{name}.index"

    sam_view_cmd = ["samtools", "view", "-Sb", "-"]
    sam_sort_cmd = ["samtools", "sort", "-"]

    stderr = subprocess.DEVNULL if config.quiet else None

    for fastq_file in glob.glob(os.path.join(sequence_dir, f"*.{file_ending}")):
        base_name = os.path.basename(fastq_file)
        bam_file = os.path.join(bam_dir, f"{base_name}.bam")
        bowtie_cmd = ["bowtie2", "-x", ref, "--very-sensitive-local", "-U", fastq_file, "--qc-filter"]
        if config.input_format == "fasta":
            bowtie_cmd.append("-f")
        sam_idx_cmd = ["samtools", "index", bam_file]

        with open(bam_file, "w") as handle:
            bowtie = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE, stderr=stderr)
            sam_view = subprocess.Popen(sam_view_cmd, stdin=bowtie.stdout, stdout=subprocess.PIPE, stderr=stderr)
            sam_sort = subprocess.Popen(sam_sort_cmd, stdin=sam_view.stdout, stdout=handle, stderr=stderr)
        sam_sort.wait()
        sam_view.wait()
        bowtie.wait()

        if bowtie.returncode != 0 or sam_view.returncode != 0 or sam_sort.returncode != 0:
            config._failed.add(bam_file)
            continue

        subprocess.check_call(sam_idx_cmd, stderr=stderr)


def _extract_if_zip(tmpdir: str, config: CSCConfig) -> str:
    """Extract the reads from a zipfile if input is indeed a zip file."""
    if os.path.isdir(config.reads):
        return config.reads
    else:
        extracted_dir = os.path.join(tmpdir, f"{config.input_format}s")
        os.makedirs(extracted_dir)
        with zipfile.ZipFile(config.reads) as zip_file:
            files = [finfo for finfo in zip_file.infolist() if finfo.filename.endswith(f".{config.input_format}")]
            for extract_file in files:
                zip_file.extract(extract_file, extracted_dir)
        return extracted_dir


def check_variants(tmpdir, config):
    bam_dir = os.path.join(tmpdir, "bams")
    outfile = open(os.path.join(config.outdir, "results.csv"), "w")

    variants = REGIONS.keys()
    columns = ["sample"]
    columns.extend(variants)
    columns.append("comment")

    print(*columns, sep=",", file=outfile)
    for bam_file in sorted(glob.glob(os.path.join(bam_dir, "*.bam"))):
        base_name = os.path.basename(bam_file)
        sample_id = base_name.split(".")[0]
        parts = [sample_id]
        found_mutations = set()

        if bam_file in config._failed:
            for variant in variants:
                parts.append("NA")
            parts.append("read failed to align")
            print(*parts, sep=",", file=outfile)
            continue

        for variant in variants:
            region = REGIONS[variant]
            try:
                before, after, quality = call_variant(config.reference, bam_file, region)
                if before == after:
                    parts.append("0")
                elif after == variant[-1]:
                    parts.append("1")
                    found_mutations.add(variant)
                else:
                    if config.show_unexpected:
                        parts.append(f"{before}{variant[1:-1]}{after}")
                    else:
                        parts.append("0")
            except PileupFailedError:
                parts.append("NA")
            except BaseDeletedError:
                parts.append("NA")
            except:
                if config.debug:
                    shutil.copy2(bam_file, "keep")
                    print(bam_file, variant)
                raise

        comment = ""
        if "D614G" not in found_mutations:
            comment += "D614G not found; low quality sequence? "
        if config.name_variants:
            comment_parts = []
            named_variants = name_variants(found_mutations)
            for variant in named_variants.keys():
                if named_variants[variant]:
                    comment_parts.append(f"found {variant}")
                else:
                    comment_parts.append(f"possibly found {variant}")

            comment += "; ".join(comment_parts)
        if "N501Y" in found_mutations and "E484K" in found_mutations:
            comment += "; important mutations found"
        comment = comment.strip()

        parts.append(comment)

        print(*parts, sep=",", file=outfile)

    outfile.close()


def name_variants(found_mutations):
    named_variants = {}
    for known_variant, expected_changes in KNOWN_VARIANTS.items():
        expected_mutations_not_found =  expected_changes.difference(found_mutations)
        unexpected_mutations_found = found_mutations.difference(expected_changes).difference({"D614G"})
        if len(expected_mutations_not_found) == 0 and len(unexpected_mutations_found) == 0:
            named_variants[known_variant] = True
        elif len(expected_mutations_not_found) <= 1 and len(unexpected_mutations_found) <= 1:
            named_variants[known_variant] = False

    return named_variants



def call_variant(reference, bam_file, region):
    cmd = ["samtools", "mpileup", "-f", reference, "-r", region, bam_file]
    mpileup = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    result = mpileup.communicate()[0].decode("utf-8")
    before, after, quality = parse_pileup(result)

    if "*" in after:
        raise BaseDeletedError()

    before_aa = Seq(before).translate()
    after_aa = Seq(after).translate()

    return before_aa, after_aa, quality


def parse_pileup(pileup):
    lines = pileup.split("\n")
    if len(lines) < 3:
        raise PileupFailedError()
    before = ""
    after = ""
    quality = []
    for line in lines[:3]:
        parts = line.split("\t")
        if len(parts) < 6:
            raise PileupFailedError()

        before += parts[2]
        after += parts[4][0] if parts[4][0] != "." else parts[2]
        quality.append(ord(parts[5])-33)

    return before, after, quality



