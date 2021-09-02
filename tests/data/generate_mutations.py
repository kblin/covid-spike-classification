#!/usr/bin/env python3

"""Add mutations to a base sequence."""


import argparse
import pathlib
import sys

from Bio import (
    Seq,
    SeqIO,
    SeqRecord,
)


def main():
    reference_path = pathlib.Path(__file__).resolve().parent.parent.parent / 'ref' / 'NC_045512.fasta'
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", type=pathlib.Path, default=reference_path,
                        help="Reference sequence to extract amplicon from (default: %(default)s).")
    parser.add_argument("-s", "--start", type=int, default=22799,
                        help="Start position of the amplicon (1-indexed) (default: %(default)s).")
    parser.add_argument("-e", "--end", type=int, default=23847,
                        help="End position of the amplicon (1-indexed, inclusive) (default: %(default)s).")
    parser.add_argument("name",
                        help="Name to use for the new sequence.")
    parser.add_argument("mutations", nargs="*", type=Mutation.fromString,
                        help="One or more mutations in the form of <position of base to mutate>:<new base>")

    args = parser.parse_args()
    run(args.reference, args.start, args.end, args.mutations, args.name)


class Mutation:
    __slots__ = ('position', 'new_base')

    def __init__(self, position, new_base):
        self.position = position
        self.new_base = new_base

    def __str__(self):
        return f"{self.position}:{self.new_base}"

    def __repr__(self):
        return f"Mut({self.position}:{self.new_base})"

    @classmethod
    def fromString(cls, string_input):
        if ":" not in string_input:
            raise ValueError(f"Invalid string format: {string_input!r}")
        pos_str, new_base = string_input.split(":", 1)
        position = int(pos_str) - 1
        if new_base not in ('A', 'C', 'G', 'T'):
            raise ValueError(f"Invalid new base {new_base!r}, must be one of 'ATGC'")

        return cls(position, new_base)


def run(reference, start, end, mutations, name):
    ref = SeqIO.read(reference, "fasta")
    seq = str(ref.seq)
    for mutation in mutations:
        seq = mutate(seq, mutation)

    amplicon = SeqRecord.SeqRecord(Seq.Seq(seq[start - 1:end]), id=name, name=name, description="")
    SeqIO.write([amplicon], sys.stdout, "fasta")


def mutate(seq: str, mutation: Mutation) -> str:
    """Mutate the base at a position given position in seq to another base"""
    before = seq[:mutation.position]
    after = seq[mutation.position + 1:]

    return before + mutation.new_base + after


if __name__ == "__main__":
    main()
