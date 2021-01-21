#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile

from .config import CSCConfig

from .core import (
    REGIONS,
    basecall,
    map_reads,
    check_variants
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-D", "--datadir", default=os.path.join(os.getcwd(), "data"),
                        help="Directory containing the ab1 files to call variants on (default: %(default)s).")
    parser.add_argument("-r", "--reference", default=os.path.join(os.getcwd(), "ref", "NC_045512.fasta"),
                        help="Reference FASTA file to use (default: %(default)s).")
    parser.add_argument("-o", "--outfile", type=argparse.FileType("w", encoding="utf-8"),
                        default=sys.stdout,
                        help="File to write result CSV to (default: stdout).")
    parser.add_argument("-q", "--quiet", action="store_true", default=False,
                        help="Suppress noisy output from the tools run")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Debug mode: Keep bam file around when the parsing crashes")
    parser.add_argument("--show-unexpected", action="store_true", default=False,
                        help="Show unexpected mutations instead of reporting 'no known mutation'")
    args = parser.parse_args()

    config = CSCConfig.from_args(args)

    with tempfile.TemporaryDirectory() as tmpdir:
        basecall(tmpdir, config)
        map_reads(tmpdir, config)
        print("sample", *REGIONS.keys(), sep=",", file=args.outfile)
        check_variants(tmpdir, config)


if __name__ == "__main__":
    main()
