"""Core functions for covid spike classification."""

import glob
import os
import shutil
import subprocess
import zipfile
from js import BioLib


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
    BioLib.set_execution_progress(5)
    fastq_dir = os.path.join(tmpdir, "fastqs")
    os.makedirs(fastq_dir)

    if os.path.isdir(config.reads):
        ab1_dir = config.reads
    else:
        ab1_dir = os.path.join(tmpdir, "ab1s")
        os.makedirs(ab1_dir)
        with zipfile.ZipFile(config.reads) as sanger_zip:
            ab1_files = [finfo for finfo in sanger_zip.infolist() if finfo.filename.endswith(".ab1")]
            for sanger_file in ab1_files:
                sanger_zip.extract(sanger_file, ab1_dir)

    sanger_files = glob.glob(os.path.join(ab1_dir, "*.ab1"))
    sanger_files_len = len(sanger_files)
    for idx, sanger_file in enumerate(sanger_files):
        BioLib.set_execution_progress(5 + (10 * (idx / sanger_files_len)))
        base_name = os.path.basename(sanger_file)
        fastq_file = os.path.join(fastq_dir, f"{base_name}.fastq")
        print("", file=open(fastq_file, 'w'))

        cmd = ["tracy", "basecall", "-f", "fastq", "-o", fastq_file, sanger_file]
        kwargs = {}
        if config.quiet:
            kwargs["stdout"] = subprocess.DEVNULL
            kwargs["stderr"] = subprocess.DEVNULL
        result = BioLib.call_task(cmd[0], b'', cmd[1:], [sanger_file, '/wasm/tracy.wasm', fastq_file], True)

    # copy all content of fastq_dir to config.outdir
    for path in os.listdir(fastq_dir):
        source_path = os.path.join(fastq_dir, path)
        destination_path = os.path.join(config.outdir, path)
        if os.path.isdir(source_path):
            shutil.copytree(source_path, destination_path)
        else:
            shutil.copy2(source_path, destination_path)


def map_reads(tmpdir, config):
    fastq_dir = os.path.join(tmpdir, "fastqs")
    bam_dir = os.path.join(tmpdir, "bams")
    os.makedirs(bam_dir)

    # ditch the .fasta file ending
    name, _ = os.path.splitext(config.reference)
    ref = f"{name}.index"
    stderr = subprocess.DEVNULL if config.quiet else None

    fqs = glob.glob(os.path.join(fastq_dir, "*.fastq"))
    fqs_length = len(fqs)

    for idx, fastq_file in enumerate(fqs):
        BioLib.set_execution_progress(15 + (25 * (idx / fqs_length)))
        base_name = os.path.basename(fastq_file)
        bam_file = os.path.join(bam_dir, f"{base_name}.bam")
        bai_file = os.path.join(bam_dir, f"{base_name}.bam.bai")
        sam_view_result_file = '/sam_view_result.bam'

        # Create output files on disk
        print("", file=open(sam_view_result_file, 'w'))
        print("", file=open(bam_file, 'w'))
        print("", file=open(bai_file, 'w'))

        sam_view_cmd = ["samtools", "view", "-Sb", "-o", "sam_view_result.bam", "-"]
        sam_sort_cmd = ["samtools", "sort", "-o", bam_file, "sam_view_result.bam"]
        bowtie_cmd = ["bowtie2", "-x", ref, "--very-sensitive-local", "-U", fastq_file, "--qc-filter"]
        sam_idx_cmd = ["samtools", "index", bam_file]


        bowtie = BioLib.call_task(bowtie_cmd[0], "", bowtie_cmd[1:],
                                  ["/ref", fastq_file, '/wasm/bowtie2-align-s.wasm'], True)
        sam_view = BioLib.call_task(sam_view_cmd[0], bytes(bowtie.stdout), sam_view_cmd[1:],
                                    [sam_view_result_file, '/wasm/samtools.wasm'], True)
        sam_sort = BioLib.call_task(sam_sort_cmd[0], "", sam_sort_cmd[1:],
                                    [sam_view_result_file, bam_file, "/wasm/samtools.wasm"], True)
        result = BioLib.call_task(sam_idx_cmd[0], "", sam_idx_cmd[1:],
                                  [bam_file, "/wasm/samtools.wasm", bai_file], True)


def check_variants(tmpdir, config):
    bam_dir = os.path.join(tmpdir, "bams")
    outfile = open(os.path.join(config.outdir, "results.csv"), "w")

    variants = REGIONS.keys()
    columns = ["sample"]
    columns.extend(variants)
    columns.append("comment")

    print(*columns, sep=",", file=outfile)
    bams = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))
    bams_length = len(bams)
    for idx, bam_file in enumerate(bams):
        BioLib.set_execution_progress(40 + (60 * (idx/bams_length)))
        base_name = os.path.basename(bam_file)
        sample_id = base_name.split(".")[0]
        parts = [sample_id]
        found_mutations = set()

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

            comment += ";".join(comment_parts)
        comment = comment.strip()
        if not comment:
            comment = "NA"

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
    result = BioLib.call_task(cmd[0], b"", cmd[1:],
                              [bam_file, reference, bam_file+'.bai', "/wasm/samtools.wasm"], True)
    before, after, quality = parse_pileup(bytes(result.stdout).decode())

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
