#!/usr/bin/env python3

import argparse
import datetime
import os
import shutil
import tempfile

from .config import CSCConfig

from .core import (
    basecall,
    map_reads,
    check_variants
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reads",
                        help="A zip file or directory containing the ab1 files to call variants on.")
    parser.add_argument("-r", "--reference", default=os.path.join(os.getcwd(), "ref", "NC_045512.fasta"),
                        help="Reference FASTA file to use (default: %(default)s).")
    parser.add_argument("-i", "--input-format", choices=["ab1", "fasta", "fastq"], default="ab1",
                        help="Select which input format to expect. Choices: %(choices)s. default: %(default)s")
    parser.add_argument("-o", "--outdir",
                        default=datetime.datetime.now().strftime("%Y-%m-%d"),
                        help="File to write result CSV and fastq files to (default: %(default)s).")
    parser.add_argument("-q", "--quiet", action="store_true", default=False,
                        help="Suppress noisy output from the tools run")
    parser.add_argument("-s", "--stdout", action="store_true", default=False,
                        help="Print results to stdout in addition to writing them to disk")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Debug mode: Keep bam file around when the parsing crashes")
    parser.add_argument("--show-unexpected", action="store_true", default=False,
                        help="Show unexpected mutations instead of reporting 'no known mutation'")
    parser.add_argument("--silence-warnings", action="store_true", default=False,
                        help="Silence D614G warnings.")
    parser.add_argument("-z", "--zip-results", action="store_true", default=False,
                        help="Create a zipfile from the output directory instead of the output directory.")
    args = parser.parse_args()

    config = CSCConfig.from_args(args)

    os.makedirs(args.outdir, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        basecall(tmpdir, config)
        map_reads(tmpdir, config)
        check_variants(tmpdir, config)
        if config.stdout:
            outfile = open(os.path.join(config.outdir, "results.csv"), "r")
            print(outfile.read())
        if config.zip_results:
            shutil.make_archive(config.outdir, "zip", root_dir=config.outdir)
            shutil.rmtree(config.outdir, ignore_errors=True)


if __name__ == "__main__":
    main()
