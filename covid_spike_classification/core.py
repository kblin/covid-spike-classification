"""Core functions for covid spike classification."""

import glob
import os
import shutil
import subprocess
import sys
import zipfile

from Bio.Seq import Seq

REGIONS = {
    "N439K": "NC_045512:22877-22879",
    "Y453F": "NC_045512:22919-22921",
    "E484K": "NC_045512:23012-23014",
    "N501Y": "NC_045512:23063-23065",
    "P681H": "NC_045512:23603-23605",
    "D614G": "NC_045512:23402-23404",
    "H655Y": "NC_045512:23525-23527",
}


class PileupFailedError(RuntimeError):
    pass


class BaseDeletedError(RuntimeError):
    pass

def basecall(tmpdir, config):
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

    for sanger_file in glob.glob(os.path.join(ab1_dir, "*.ab1")):
        base_name = os.path.basename(sanger_file)
        cmd = ["tracy", "basecall", "-f", "fastq", "-o", os.path.join(fastq_dir, f"{base_name}.fastq"), sanger_file]
        kwargs = {}
        if config.quiet:
            kwargs["stdout"] = subprocess.DEVNULL
            kwargs["stderr"] = subprocess.DEVNULL
        subprocess.check_call(cmd, **kwargs)

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

    sam_view_cmd = ["samtools", "view", "-Sb", "-"]
    sam_sort_cmd = ["samtools", "sort", "-"]

    stderr = subprocess.DEVNULL if config.quiet else None

    for fastq_file in glob.glob(os.path.join(fastq_dir, "*.fastq")):
        base_name = os.path.basename(fastq_file)
        bam_file = os.path.join(bam_dir, f"{base_name}.bam")
        bowtie_cmd = ["bowtie2", "-x", ref, "--very-sensitive-local", "-U", fastq_file, "--qc-filter"]
        sam_idx_cmd = ["samtools", "index", bam_file]

        with open(bam_file, "w") as handle:
            bowtie = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE, stderr=stderr)
            sam_view = subprocess.Popen(sam_view_cmd, stdin=bowtie.stdout, stdout=subprocess.PIPE, stderr=stderr)
            sam_sort = subprocess.Popen(sam_sort_cmd, stdin=sam_view.stdout, stdout=handle, stderr=stderr)
            sam_sort.wait()

        subprocess.check_call(sam_idx_cmd, stderr=stderr)


def check_variants(tmpdir, config):
    bam_dir = os.path.join(tmpdir, "bams")
    outfile = open(os.path.join(config.outdir, "results.csv"), "w")

    print("sample", *REGIONS.keys(), sep=",", file=outfile)
    for bam_file in sorted(glob.glob(os.path.join(bam_dir, "*.bam"))):
        base_name = os.path.basename(bam_file)
        sample_id = base_name.split(".")[0]
        parts = [sample_id]
        for variant, region in REGIONS.items():
            try:
                before, after, quality = call_variant(config.reference, bam_file, region)
                if before == after:
                    parts.append("n")
                elif after == variant[-1]:
                    parts.append("y")
                else:
                    if config.show_unexpected:
                        parts.append(f"{before}{variant[1:-1]}{after}")
                    else:
                        parts.append("n")
            except PileupFailedError:
                parts.append("NA")
            except BaseDeletedError:
                parts.append("0")
            except:
                if config.debug:
                    shutil.copy2(bam_file, "keep")
                    print(bam_file, variant)
                raise
        print(*parts, sep=",", file=outfile)

    outfile.close()


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



