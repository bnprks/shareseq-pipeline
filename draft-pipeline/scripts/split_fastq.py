#!/usr/bin/env python3

# Authors: Ben Parks
# Last updated: 9/26/22

# Split a fastq file into multiple files with a fixed number of reads per file
# Outputs will be named 1.fastq.zst, 2.fastq.zst ... in the specified directory
# Here zstandard (zst) compression is used on outputs for much faster read + write

import argparse
import gzip
import os
import subprocess

from typing import *
from typing.io import *

def main():
    parser = argparse.ArgumentParser(description="Split a fastq file into multiple fixed-length chunks")
    parser.add_argument("input_fastq", type=str, help="Fastq file to be split")
    parser.add_argument("output_path", type=str, help="Directory to put output chunks")
    parser.add_argument("--reads", type=int, help="Number of reads per chunk")

    args = parser.parse_args()

    input_file = gzip.open(args.input_fastq)
    
    chunk_number = 0
    lines_written = args.reads * 4
    os.mkdir(args.output_path)
    output_file = None
    
    for line in input_file:
        if lines_written >= args.reads*4:
            lines_written = 0
            if output_file:
                output_file.close()
            chunk_number = chunk_number + 1
            output_file = open_zstd(os.path.join(args.output_path, f"{chunk_number}.fastq.zst"))
        output_file.write(line)
        lines_written += 1

def open_zstd(path: str) -> BinaryIO:
    ps = subprocess.Popen(('zstd', '-c', '--fast=1'), stdout=open(path, 'wb'), stdin=subprocess.PIPE)
    return ps.stdin

if __name__ == "__main__":
    main()
