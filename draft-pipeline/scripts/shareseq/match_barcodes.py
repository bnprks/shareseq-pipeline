#!/usr/bin/env python3

# Authors: Ben Parks
# Last updated: 9/26/22

# Match barcodes from SHARE-seq data
# Input: 
#   - Uncompressed fastqs for R1, R2, where I1 + I2 sequences are appended to read names as format :[I1]+[I2]
#   - TSV of valid barcodes and labels
#   - List of I2 sequences to output as ATAC library
#   - List of I2 sequences to output as RNA library
# Output: 
#   - Uncompressed compressed fastqs for R1 and R2 split between ATAC and RNA.
#     Read name has the cell barcode ID appended rather than raw sequence
#   - QC file in json format, containing total reads and number of reads matching
#     each barcode with 0 or 1 mismatches

# Details: 
# - Performs up to 1bp correction in each barcode
# - 3 Rounds of barcodes are matched against the I1, and output in the read name
# - The I2 sequence is only used for separating ATAC + RNA libraries, and is not
#   incorporated into the read name.

import argparse
import collections
import itertools
import json

from typing import *
from typing.io import *

def main():
    parser = argparse.ArgumentParser(description="Match SHARE-seq barcodes")
    parser.add_argument("--R1_in", type=str, help="Fastq R1 input")
    parser.add_argument("--R2_in", type=str, help="Fastq R2 input")

    parser.add_argument("--R1_out", required=True, type=str, help="Fastq R1 output")
    parser.add_argument("--R2_out", required=True, type=str, help="Fastq R2 output")

    parser.add_argument("--BC1", required=True, type=str, help="TSV of BC1 sequences")
    parser.add_argument("--BC2", required=True, type=str, help="TSV of BC2 sequences")
    parser.add_argument("--BC3", required=True, type=str, help="TSV of BC3 sequences")

    parser.add_argument("--BC1_pos", default=15, type=str, help="Position of BC1 (0-based)")
    parser.add_argument("--BC2_pos", default=53, type=str, help="Position of BC2 (0-based)")
    parser.add_argument("--BC3_pos", default=91, type=str, help="Position of BC3 (0-based)")

    parser.add_argument("--json_stats", required=True, type=str, help="Path for json matching stats output")

    args = parser.parse_args()
    
    # Maximum number of total mismatches (note that each barcode can only have up to 1 mismatch)
    max_mismatches = 3
    
    # Prep the I1 index sequence lookups
    barcodes = []
    positions = []
    bc_lengths = []
    for i in [1,2,3]:
        seqs = read_bc(vars(args)[f"BC{i}"])
        seqs = add_mismatches(seqs)
        barcodes.append(seqs)
        positions.append(vars(args)[f"BC{i}_pos"])
        bc_lengths.append(len(next(iter(seqs))))

    # Prep stats counters
    bc_match_stats = [collections.defaultdict(int) for i in range(3)]
    mismatch_counts = [0 for i in range(max_mismatches + 2)]

    # Prep the output files
    R1_out = open(args.R1_out, "wb")
    R2_out = open(args.R2_out, "wb")

    json_output = open(args.json_stats, "w")

    # Loop through reads and perform matching
    current_match = ["", "", ""]
    for r1, r2 in zip(read_fastq(args.R1_in), read_fastq(args.R2_in)):
        bc1_pos = r1.name.rfind(b":") + 1
        
        # For each of the 3 I1 barcodes:
        # Check if the barcode at that position matches
        total_mismatches = 0
        for i in range(3):
            offset = bc1_pos + positions[i]
            bc = r1.name[offset : offset+bc_lengths[i]]
            if bc not in barcodes[i]:
                total_mismatches = 4
            else:
                name, mismatches = barcodes[i][bc]
                total_mismatches += mismatches
                bc_match_stats[i][(name, mismatches)] += 1
                current_match[i] = name
        
        total_mismatches = min(total_mismatches, max_mismatches + 1)
        mismatch_counts[total_mismatches] += 1
        
        if total_mismatches <= max_mismatches:
            write_fastq(R1_out, r1, current_match)
            write_fastq(R2_out, r2, current_match)

    # Output stats
    stats = {}
    for i in [1,2,3]:
        stats[f"BC{i}"] = format_stats(bc_match_stats[i-1])
    stats["total_mismatch_histogram"] = mismatch_counts
    json.dump(stats, json_output)


Read = collections.namedtuple("Read", ["name", "seq", "qual"])
def read_fastq(file: str) -> Iterator[Read]:
    """Yield tuples of the read from a fastq file iterator"""
    file_iter = open(file, "rb")
    while True:
        name = next(file_iter, None)
        if name is None:
            return None
        seq = next(file_iter, None)
        _ = next(file_iter, None)
        qual = next(file_iter, None)
        yield Read(name, seq, qual)

def write_fastq(file: BinaryIO, read: Read, barcodes: List[bytes]):
    # Write name with barcode appended
    file.write(read.name[:read.name.find(b" ")])
    file.write(b" CB:Z:")
    for b in barcodes[:-1]:
        file.write(b)
        file.write(b"_")
    file.write(barcodes[-1])
    file.write(b"\n")

    file.write(read.seq)
    file.write(b"+\n")
    file.write(read.qual)


def read_bc(tsv_path: str) -> Dict[bytes, bytes]:
    """Return dictionary from sequence -> label for barcodes in a tsv format"""
    lines = open(tsv_path, "rb").readlines()
    assert lines[0] == b"Name\tSequence\n"
    seqs = {}
    for l in lines[1:]:
        name, seq = l.strip().split(b"\t")
        seq = seq.upper()
        seqs[seq] = name
    assert len(set(len(seq) for seq in seqs.keys())) == 1
    return seqs


def single_mismatches(seq: bytes) -> Iterator[bytes]:
    """Iterate through 1bp mismatches of a sequence"""
    for idx, base in itertools.product(range(len(seq)), [b"A", b"T", b"G", b"C", b"N"]):
        if base == seq[idx:idx+1]:
            continue
        yield seq[:idx] + base + seq[idx+1:]


def add_mismatches(seqs: Dict[bytes, bytes]) -> Dict[bytes, Tuple[bytes, int]]:
    """
    Input: dictionary from sequence -> label
    Output: dictionary from sequence -> (label, mismatches)
    Where all 1bp mismatches have been added to the dictionary
    """
    output = {}
    for seq, name in seqs.items():
        assert seq not in output
        output[seq] = (name, 0)
        for mutant in single_mismatches(seq):
            assert mutant not in seqs
            if mutant in output:
                output[mutant] = None
            else:
                output[mutant] = (name, 1)
    return {k:v for k,v in output.items() if v is not None}


def format_stats(stats: Dict[Tuple[bytes, int], int]) -> Dict[str, Dict[str, int]]:
    """Input: stats dictionary with keys (name, mismatches) 
    Output: Dictionary of "exact_match" and "1bp_mismatch" counts, 
    where each entry is a dictionary from barcode to count
    """
    exact_matches = {name.decode(): v for (name, mismatches), v in stats.items() if mismatches == 0}
    mismatches = {name.decode(): v for (name, mismatches), v in stats.items() if mismatches == 1}
    return {
        "exact_match": exact_matches,
        "1bp_mismatch": mismatches,
    }


if __name__ == "__main__":
    main()
