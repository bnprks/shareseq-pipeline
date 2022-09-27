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
    parser = argparse.ArgumentParser(description="Split a fastq file into multiple fixed-length chunks")
    parser.add_argument("--R1_in", type=str, help="Fastq R1 input")
    parser.add_argument("--R2_in", type=str, help="Fastq R2 input")

    parser.add_argument("--R1_ATAC", required=True, type=str, help="Fastq R1 ATAC output")
    parser.add_argument("--R2_ATAC", required=True, type=str, help="Fastq R2 ATAC output")

    parser.add_argument("--R1_RNA", required=True, type=str, help="Fastq R1 RNA output")
    parser.add_argument("--R2_RNA", required=True, type=str, help="Fastq R2 RNA output")

    parser.add_argument("--ATAC_I2", required=True, type=str, nargs='+', help="I2 barcode sequences for ATAC libraries")
    parser.add_argument("--RNA_I2", required=True, type=str, nargs='+', help="I2 barcode sequences for RNA libraries")

    parser.add_argument("--BC1", required=True, type=str, help="TSV of BC1 sequences")
    parser.add_argument("--BC2", required=True, type=str, help="TSV of BC2 sequences")
    parser.add_argument("--BC3", required=True, type=str, help="TSV of BC3 sequences")

    parser.add_argument("--BC1_pos", default=15, type=str, help="Position of BC1 (0-based)")
    parser.add_argument("--BC2_pos", default=53, type=str, help="Position of BC2 (0-based)")
    parser.add_argument("--BC3_pos", default=91, type=str, help="Position of BC3 (0-based)")

    parser.add_argument("--json_stats", required=True, type=str, help="Path for json matching stats output")

    args = parser.parse_args()
    
    # Prep the I2 index sequence lookup
    i2_length = len(args.ATAC_I2[0])
    assert all(len(bc) == i2_length for bc in args.ATAC_I2 + args.RNA_I2)

    i2_atac_seqs = {bc.encode(): bc.encode() for bc in args.ATAC_I2}
    i2_atac_seqs = add_mismatches(i2_atac_seqs)

    i2_rna_seqs = {bc.encode(): bc.encode() for bc in args.RNA_I2}
    i2_rna_seqs = add_mismatches(i2_rna_seqs)
    
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

    # Prep stats dictionaries
    bc_match_stats = [collections.defaultdict(int) for i in range(3)]
    atac_match_stats = collections.defaultdict(int)
    rna_match_stats = collections.defaultdict(int)

    # Prep the output files
    ATAC_R1 = open(args.R1_ATAC, "wb")
    ATAC_R2 = open(args.R2_ATAC, "wb")
    RNA_R1 = open(args.R1_RNA, "wb")
    RNA_R2 = open(args.R2_RNA, "wb")

    json_output = open(args.json_stats, "w")

    # Loop through reads and perform matching
    current_match = ["", "", ""]
    for r1, r2 in zip(read_fastq(args.R1_in), read_fastq(args.R2_in)):
        bc1_pos = r1.name.rfind(b":") + 1
        bc2_pos = r1.name.rfind(b"+") + 1
        
        # For each of the 3 I1 barcodes:
        # Check if the barcode at that position matches
        has_match = True
        for i in range(3):
            offset = bc1_pos + positions[i]
            bc = r1.name[offset : offset+bc_lengths[i]]
            if bc not in barcodes[i]:
                has_match = False
            else:
                name, mismatches = barcodes[i][bc]
                bc_match_stats[i][(name, mismatches)] += 1
                current_match[i] = name
        
        # Check the I2 barcode and output
        bc2 = r1.name[bc2_pos : bc2_pos+i2_length]
        if bc2 in i2_atac_seqs:
            atac_match_stats[i2_atac_seqs[bc2]] += 1
            if has_match:
                write_fastq(ATAC_R1, r1, current_match)
                write_fastq(ATAC_R2, r2, current_match)

        if bc2 in i2_rna_seqs:
            rna_match_stats[i2_rna_seqs[bc2]] += 1
            if has_match:
                write_fastq(RNA_R1, r1, current_match)
                write_fastq(RNA_R2, r2, current_match)

    # Output stats
    stats = {}
    for i in [1,2,3]:
        stats[f"BC{i}"] = format_stats(bc_match_stats[i-1])
    stats["ATAC_BC"] = format_stats(atac_match_stats)
    stats["RNA_BC"] = format_stats(rna_match_stats)
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
    file.write(b" ")
    for b in barcodes[:-1]:
        file.write(b)
        file.write(b"_")
    file.write(barcodes[-1])
    file.write(b"\n")

    file.write(read.seq)
    file.write(b"+\n")
    file.write(read.qual)


complement = bytes.maketrans(b"ATGC", b"TACG")
def reverse_complement(seq):
    return seq.translate(complement)[::-1]


def read_bc(tsv_path: str) -> Dict[bytes, bytes]:
    """Return dictionary from sequence -> label for barcodes in a tsv format"""
    lines = open(tsv_path, "rb").readlines()
    assert lines[0] == b"Name\tSequence\n"
    seqs = {}
    for l in lines[1:]:
        name, seq = l.strip().split(b"\t")
        seq = reverse_complement(seq.upper())
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
