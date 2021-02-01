#!/usr/bin/env python3

import argparse
import datetime
import os
import shutil
import tempfile
from js import BioLib

from covid_spike_classification.config import CSCConfig

from covid_spike_classification.core import (
    REGIONS,
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
    parser.add_argument("-n", "--name-variants", action="store_true", default=False,
                        help="Add a column naming known variants")
    parser.add_argument("-z", "--zip-results", action="store_true", default=False,
                        help="Create a zipfile from the output directory instead of the output directory.")
    args = parser.parse_args()

    config = CSCConfig.from_args(args)

    os.makedirs(args.outdir, exist_ok=True)

    # Build indices
    BioLib.call_task(
        "bowtie2-build",
        b'',
        "/ref/NC_045512.fasta /ref/NC_045512.index",
        ["/ref", "/wasm/bowtie2-build-s.wasm"],
        True)

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
