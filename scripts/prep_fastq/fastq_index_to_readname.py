
# Add I1 + I2 sequences to the end of R1 + R2 sequences
# Output read names will be like:
# @A00509:608:HNT5JDSX3:1:2154:1063:1016 1:N:0:NCGGACGATCATG+TTTCNAGC
# Where the last sequence bits are I1 then I2

# Authors: Ben Parks
# Last updated: 9/27/22

# Usage: cat input.tsv | python dedup_fragments.py > output.tsv

import argparse
import collections
import os
import subprocess

from typing import *
from typing.io import *

def main():
    parser = argparse.ArgumentParser(description="Add I1+I2 sequences to the end of R1 + R2")
    parser.add_argument("--R1", type=str, help="R1 input")
    parser.add_argument("--R2", type=str, help="R2 input")
    parser.add_argument("--I1", type=str, help="I1 input")
    parser.add_argument("--I2", type=str, help="I2 input")
    parser.add_argument("--R1_out", type=str, help="R1 output")
    parser.add_argument("--R2_out", type=str, help="R2 output")
    args = parser.parse_args()

    R1_out = write_gzip(args.R1_out, level=4)
    R2_out = write_gzip(args.R2_out, level=4)

    inputs = zip(read_fastq(args.R1), read_fastq(args.R2), read_fastq(args.I1), read_fastq(args.I2))
    for r1, r2, i1, i2 in inputs:
        write_fastq(R1_out, r1, b"1:N:0:", i1.seq[:-1], i2.seq[:-1])
        write_fastq(R2_out, r2, b"2:N:0:", i1.seq[:-1], i2.seq[:-1])


def read_gzip(filename: str) -> BinaryIO:
    """Read a file, using an external process to decompress gzip if needed"""
    if os.path.isfile(filename) and not filename.endswith(".gz"):
        return open(filename, 'rb')
    else:
        if not os.path.isfile(filename):
            raise RuntimeError(f"File does not exist: {filename}")
        ps = subprocess.Popen(('gzip', '-d', '-c', filename), stdout=subprocess.PIPE)
        return ps.stdout


def write_gzip(filename: str, level: int = 4) -> BinaryIO:
    """Write a file, using an external process to compress gzip if needed"""
    if filename[-3:] == ".gz":
        ps = subprocess.Popen(('gzip', '-c', f'-{level}'), stdout=open(filename, 'wb'), stdin=subprocess.PIPE)
        return ps.stdin
    else:
        return open(filename, 'wb')

Read = collections.namedtuple("Read", ["name", "seq", "qual"])
def read_fastq(file: str) -> Iterator[Read]:
    """Yield tuples of the read from a fastq file iterator"""
    file_iter = read_gzip(file)
    while True:
        name = next(file_iter, None)
        if name is None:
            return None
        seq = next(file_iter, None)
        _ = next(file_iter, None)
        qual = next(file_iter, None)
        yield Read(name, seq, qual)

def write_fastq(file, read: Read, prefix, i1_seq, i2_seq):
    # Write name with barcode appended
    file.write(read.name[:read.name.find(b" ")+1])
    file.write(prefix)
    file.write(i1_seq)
    file.write(b"+")
    file.write(i2_seq)
    file.write(b"\n")

    file.write(read.seq)
    file.write(b"+\n")
    file.write(read.qual)

if __name__ == "__main__":
    main()
