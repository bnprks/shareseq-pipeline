#!/usr/bin/env python3

# Authors: Ben Parks
# Last updated: 9/26/22

# Run bcl2fastq for a share-seq library
# Inputs:
#   - bcl run directory
#   - snakemake config object with sublibrary barcodes
# Outputs:
#   - I1, I2, R1, R2 fastq.gz files for ATAC and RNA in their respective directories

# Details:
# Because bcl2fastq doesn't like to demultiplex on only some of the index reads,
# we use a bases mask to pretend that I1 is not an index read.
# After that, we rename files so that we keep the normall I1/I2/R1/R2 notation
from pathlib import Path
import subprocess

import re


def main():
    out_dir = snakemake.output.results_dir
    run_dir = snakemake.input.run_dir
    sequencing_run = snakemake.config["sequencing"][snakemake.wildcards.sequencing_run]
    assert sequencing_run["type"] == "bcl"
    
    # Make samplesheet
    with open(snakemake.output.samplesheet, "w") as out:
        out.write("[Data]\nSample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n")
        for sample_id, I2 in sequencing_run["ATAC_I2"].items():
            out.write(f"ATAC_{sample_id},,,,bc,,,{I2},,\n")
        for sample_id, I2 in sequencing_run["RNA_I2"].items():
            out.write(f"RNA_{sample_id},,,,bc,,,{I2},,\n")

    log = open(snakemake.log[0], "wb")

    # Run bcl2fastq
    command = [
        "bcl2fastq", 
        "--runfolder-dir", run_dir,
        "--sample-sheet", snakemake.output.samplesheet, 
        "--no-lane-splitting",
        "--create-fastq-for-index-reads",
        "--use-bases-mask", get_bases_mask(run_dir),
        "--output", out_dir,
        "--processing-threads", str(snakemake.threads), # Empirically, vast majority of the work is spent in processing threads, so we don't add threads for other purposes
        # "--tiles", "1_2154,2_2424,3_2304,4_1557" # Uncomment to use a subset of reads for testing
    ]
    print(" ".join(command))
    bcl2fastq = subprocess.run(command, stderr=log)
    log.close()

    if bcl2fastq.returncode != 0:
        sys.exit(bcl2fastq.returncode)

    # Rename output files to match the right file naming convention
    Path(f"{out_dir}/ATAC").mkdir()
    Path(f"{out_dir}/RNA").mkdir()
    i = 0
    read_rename = {
        "R1": "R1", 
        "R2": "I1",
        "R3": "R2",
        "I1": "I2",
    }
    for sublib_id in sequencing_run["ATAC_I2"]:
        for old_read, new_read in read_rename.items():
            bcl2fastq_path = Path(f"{out_dir}/ATAC_{sublib_id}_S{i+1}_{old_read}_001.fastq.gz")
            bcl2fastq_path.rename(f"{out_dir}/ATAC/{sublib_id}_{new_read}.fastq.gz")
        i += 1
    for sublib_id in sequencing_run["RNA_I2"]:
        for old_read, new_read in read_rename.items():
            bcl2fastq_path = Path(f"{out_dir}/RNA_{sublib_id}_S{i+1}_{old_read}_001.fastq.gz")
            bcl2fastq_path.rename(f"{out_dir}/RNA/{sublib_id}_{new_read}.fastq.gz")
        i += 1

def get_bases_mask(run_dir):
    """Parse out a bases mask from the run_dir, switching I1 to be treated like a read"""
    # Do some regular expression hunting to parse the read cycle counts
    raw = [
        l.strip() 
        for l in open(Path(run_dir) / "RunInfo.xml")
        if "<Read Number=" in l
    ]
    assert len(raw) == 4
    parsed = [
        re.match('<Read Number="(\d)" NumCycles="(\d+)" IsIndexedRead="(Y|N)"/>', read)
        for read in raw
    ]
    # Carefully check our parsing
    for i, match in enumerate(parsed):
        assert match is not None
        assert int(match.group(1)) == i+1
        if i in [0,3]:
            assert match.group(3) == "N"
        else:
            assert match.group(3) == "Y"
    # Make the bases mask
    return (
        f"Y{parsed[0].group(2)},"
        f"Y{parsed[1].group(2)},"
        f"I{parsed[2].group(2)},"
        f"Y{parsed[3].group(2)}"
    )

if __name__ == "__main__":
    main()