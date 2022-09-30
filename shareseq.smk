# Main shareseq pipeliine
# Author: Ben Parks
# Last Modified: 9/28/22

# Primary outputs:
# - ATAC/fragments.tsv.gz: 10x-compatible fragment file (though end coordinates are +1bp relative to 10x)

import os

workdir: config["output_dir"]

chunk_size = config["chunk_size"]

# Fix up any config keys that are not strings
def string_only_keys(data):
    """Recursively edit a dictionary to convert all keys to string types"""
    assert isinstance(data, dict)
    for k, v in list(data.items()):
        if not isinstance(k, str):
            data[str(k)] = v
            del data[k]
        if isinstance(v, dict):
            string_only_keys(v)
string_only_keys(config)

def get_sequencing_paths(assay):
    """Get a list of all paths for sublibraries for the current sample"""
    assert assay in ["ATAC", "RNA"]
    sequencing_paths = []
    for run_id, run in config["sequencing"].items():
        assert run["type"] in ["bcl", "fastq"]
        if run["type"] == "bcl":
            sequencing_paths += [
                f"{assay}/{run_id}/{sublib_id}" for sublib_id in run[f"{assay}_I2"]
            ]
        elif run["type"] == "fastq":
            sequencing_paths.append(f"{assay}/{run_id}")
        else:
            assert False
    return sequencing_paths

# Confirm that we have read counts for all the input sequences
for sequencing_path in get_sequencing_paths("ATAC") + get_sequencing_paths("RNA"):
    if not os.path.exists(f"{sequencing_path}/read_count.txt"):
        raise RuntimeError(f"Must run prep_fastq.smk; missing read counts for: {sequencing_path}")
del sequencing_path

wildcard_constraints:
    chunk = "\d+", # Chunk is a number
    sequencing_path = "(ATAC|RNA)/([^/]+/)?[^/]+", # Sequencing path is 2-3 folders

rule all:
    input: 
        'ATAC/fragments.tsv.gz',
        "ATAC/alignment_stats.json",

#############################
### ATAC + RNA fastq processing 
#############################
def fastq_path(sequencing_path, read):
    """Take a sublibrary path and return the path to its R1 or R2 fastq"""
    assay_type, run_id = sequencing_path.split("/")[:2] 
    if config["sequencing"][run_id]["type"] == "bcl":
        return f"bcl2fastq/{sequencing_path}_{read}.fastq.zst"
    elif config["sequencing"][run_id]["type"] == "fastq":
        return config["sequencing"][run_id][assay_type.upper() + "_" + read]
    else:
        assert False

def fastq_decompress(sequencing_path):
    """Take a sublibrary path and return the command to decompress it"""
    assay_type, run_id = sequencing_path.split("/")[:2] 
    if config["sequencing"][run_id]["type"] == "bcl":
        return "zstd -dc"
    elif config["sequencing"][run_id]["type"] == "fastq":
        return "gzip -dc"
    else:
        assert False

# Split fastqs
rule split_fastqs:
    input:
        fastq = lambda w: fastq_path(w.sequencing_path, w.read),
        read_count = "{sequencing_path}/read_count.txt",
    output:
        chunks = temp(directory("{sequencing_path}/split_fastqs/{read}"))
    params:
        decompress = lambda w: fastq_decompress(w.sequencing_path),
        lines = chunk_size * 4,
    resources:
        runtime = 60 * 5, # Be generous on time in case of large fastqs
    threads: 3
    log: '{sequencing_path}/split_fastqs/{read}.log'
    shell: "mkdir {output.chunks} && "
          " split <({params.decompress} {input.fastq}) "
          " --numeric-suffixes=1 --lines {params.lines} "
          " --additional-suffix=.fastq.zst "
          " --filter='zstd --fast=1 -q -o $FILE' "
          " {output.chunks}/ 2> {log}"

# Perform barcode matching
rule match_barcodes:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True)
    output:
        R1 = temp("{sequencing_path}/{chunk}/01_match_barcodes_R1.fastq.zst"),
        R2 = temp("{sequencing_path}/{chunk}/01_match_barcodes_R2.fastq.zst"),
        stats = "{sequencing_path}/{chunk}/qc_stats/01_match_barcodes.json",
    params:
        script = srcdir("scripts/shareseq/match_barcodes.py"),
        R1_in = "{sequencing_path}/split_fastqs/R1/{chunk}.fastq.zst",
        R2_in = "{sequencing_path}/split_fastqs/R2/{chunk}.fastq.zst",
        BC1 = srcdir("config/barcodes/Round1.tsv"),
        BC2 = srcdir("config/barcodes/Round2.tsv"),
        BC3 = srcdir("config/barcodes/Round3.tsv"),
    threads: 2
    log: '{sequencing_path}/{chunk}/01_match_barcodes.log'
    # Here we use process substitution like >() and <() to use zstandard to 
    # compress + decompress our inputs/outputs
    shell: "python3 {params.script} "
        " --R1_in <(zstd -dc {params.R1_in}) --R2_in <(zstd -dc {params.R2_in}) "
        " --R1_out >(zstd --fast=1 -q -o {output.R1}) "
        " --R2_out >(zstd --fast=1 -q -o {output.R2}) "
        " --BC1 {params.BC1} --BC2 {params.BC2} --BC3 {params.BC3} "
        " --json_stats {output.stats} "
        " 2> {log} "

# Remove adapter ends from the raw fastq reads
rule trim_adapters:
    input: 
        R1 = rules.match_barcodes.output.R1,
        R2 = rules.match_barcodes.output.R2,
    output:
        interleaved = temp("{sequencing_path}/{chunk}/02_trim_adapters.interleaved.fastq.zst"),
        report_json = "{sequencing_path}/{chunk}/qc_stats/02_trim_adapters.json"
    threads: 4
    log: '{sequencing_path}/{chunk}/02_trim_adapters.log'
    shell: "fastp --in1 <(zstd -dc {input.R1}) --in2 <(zstd -dc {input.R2}) "
        " --adapter_sequence    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC "
        " --adapter_sequence_r2 CTGTCTCTTATACACATCTGACGCTGCCGACGA "
        " -j {output.report_json} -G -Q -L -w {threads} 2> {log} "
        " --stdout | zstd --fast=1 -q -o {output.interleaved}"

# Alternative trim command
# "SeqPurge -a1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a2 CTGTCTCTTATACACATCTGACGCTGCCGACGA "
#     " -qcut 0 -ncut 0 "
#     " -threads {threads} -out1 {output.R1} -out2 {output.R2} "
#     " -in1 {input.R1} -in2 {input.R2} > {log}"

# fastp args from ENCODE sc_atac pipeline
# "fastp -i {input.fastq1_bc} -I {input.fastq2_bc} -o {output.fastq1_trim} -O {output.fastq2_trim}"
#        " -h {log.html} -j {log.json} -G -Q -L -w {threads} 2> {output.stats}"

#############################
### ATAC-specific workflow 
#############################

# Align ATAC reads with bowtie2, and filter to good quality reads only
rule atac_bowtie2:
    input: 
        fastq = rules.trim_adapters.output.interleaved
    output: 
        bam = temp('{sequencing_path}/{chunk}/03_atac_bowtie2.bam'),
    params:
        index = config["genome"]["bowtie2"]
    resources:
        runtime = min(60, 2 * config["chunk_size"] // 1_000_000) # 2 minutes-per 1M read time estimate
    threads: 16
    log: '{sequencing_path}/{chunk}/03_atac_bowtie2.log',
    shell: "bowtie2 --interleaved <(zstd -dc {input.fastq}) -x {params.index} "
           " --sam-append-comment --maxins 2000 --threads {threads} 2> {log} | "
           # -F 1804: exclude flag, exludes unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
           # -f 2: flags to require, properly aligned
           # -q 30: exlude low MAPQ, set as adjustable configuration parameter
           "samtools view -F 1804 -f 2 -q 30 -1 - > {output} "
        
# Convert bam to fragments format and sort for first pass
rule atac_convert_fragments:
    input:
        bam = rules.atac_bowtie2.output.bam
    output: 
        fragments = temp('{sequencing_path}/{chunk}/04_atac_convert_fragments.fragments.tsv'),
    params:
        script = srcdir("scripts/shareseq/bam_to_fragments.py"),
        memory = "4G",
    threads: 4
    log: '{sequencing_path}/{chunk}/04_atac_convert_fragments.log',
    shell: "python {params.script} {input} 2> {log} | "
           "LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 -t$'\\t' "
           "-S {params.memory} --parallel={threads} > {output.fragments} "

def get_chunks(sequencing_path):
    """Generate chunk IDs for a sequencing path based on read count. Adds padding 0s as needed"""
    reads = int(open(f"{sequencing_path}/read_count.txt").read())
    chunk_count = (reads + chunk_size - 1) // chunk_size
    str_len = max(2, len(str(chunk_count)))
    return [f"{i:0{str_len}d}" for i in range(1, chunk_count+1)]

rule atac_merge_chunks:
    input: lambda w: expand(rules.atac_convert_fragments.output, sequencing_path=w.sequencing_path, chunk=get_chunks(w.sequencing_path))
    output: 
        fragments = temp('{sequencing_path}/fragments.tsv'),
        compressed = '{sequencing_path}/fragments.tsv.gz',
        indexed = '{sequencing_path}/fragments.tsv.gz.tbi',
    params:
        memory = "4G",
        script = srcdir("scripts/shareseq/dedup_fragments.py")
    shell: "LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 -t$'\\t' "
        " --merge --batch-size=100 -S {params.memory} {input} | "
        " python {params.script} | "
        " tee {output.fragments} | "
        " bgzip -c > {output.compressed} && "
        " tabix --preset bed {output.compressed}"

rule atac_merge_libraries:
    input: lambda w: expand(rules.atac_merge_chunks.output.fragments, sequencing_path=get_sequencing_paths("ATAC"))
    output:
        compressed = 'ATAC/fragments.tsv.gz',
        indexed = 'ATAC/fragments.tsv.gz.tbi',
    params:
        memory = "4G",
        script = srcdir("scripts/shareseq/dedup_fragments.py")
    shell:"LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 -t$'\\t' "
        " --merge --batch-size=100 -S {params.memory} {input} | "
        " python {params.script} | "
        " bgzip -c > {output.compressed} && "
        " tabix --preset bed {output.compressed} "

rule atac_stats:
    input: rules.atac_merge_libraries.output.compressed
    output: 
        json = "ATAC/alignment_stats.json"
    params:
        script = srcdir("scripts/shareseq/collect_atac_stats.py")
    shell: "python {params.script} > {output.json}"
