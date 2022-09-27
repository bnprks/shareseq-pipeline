
# Number of reads for a chunk
chunk_size = 1_000_000

# Config structure:
# samples:
#   sample_name:
#       R1: path/to/R1.fastq.gz
#       R2: path/to/R2.fastq.gz
#       ATAC_I2: ["ATGCATA", "GCATACAA"]
#       RNA_I2:  ["ATGCATA", "CCGGATAA"]

rule all:
    input:
        expand("{sample}/fastq/match_barcodes/{assay}/{chunk}_R1.fastq.zst", 
            sample=config["samples"].keys(),
            assay=["ATAC", "RNA"],
            chunk=[1]),
        expand("{sample}/fastq/trim_adapters/{assay}/{chunk}_interleaved.fastq.zst",
            sample=config["samples"].keys(),
            assay=["ATAC", "RNA"],
            chunk=[1])

# Split fastqs
rule split_fastqs:
    input:
        fastq = lambda w: config["samples"][w.sample][w.read]
    output:
        chunks = temp(directory("{sample}/fastq/split_fastqs/{read}"))
    params:
        script = srcdir("scripts/split_fastq.py")
    log: '{sample}/logs/split_fastqs_{read}.log'
    shell: "python3 {params.script} {input} {output} --reads {chunk_size} 2> {log}"

# Perform barcode matching
rule match_barcodes:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True)
    output:
        R1_RNA = temp("{sample}/fastq/match_barcodes/RNA/{chunk}_R1.fastq.zst"),
        R2_RNA = temp("{sample}/fastq/match_barcodes/RNA/{chunk}_R2.fastq.zst"),
        R1_ATAC = temp("{sample}/fastq/match_barcodes/ATAC/{chunk}_R1.fastq.zst"),
        R2_ATAC = temp("{sample}/fastq/match_barcodes/ATAC/{chunk}_R2.fastq.zst"),
        stats = "{sample}/qc_stats/chunks/match_barcodes/{chunk}.json",
    params:
        script = srcdir("scripts/match_barcodes.py"),
        ATAC_I2 = lambda w: " ".join(config["samples"][w.sample]["ATAC_I2"]),
        RNA_I2 = lambda w: " ".join(config["samples"][w.sample]["RNA_I2"]),
        BC1 = srcdir("config/barcodes/Round1.tsv"),
        BC2 = srcdir("config/barcodes/Round2.tsv"),
        BC3 = srcdir("config/barcodes/Round3.tsv"),
        R1_in = "{sample}/fastq/split_fastqs/R1/{chunk}.fastq.zst",
        R2_in = "{sample}/fastq/split_fastqs/R2/{chunk}.fastq.zst",
    log: '{sample}/logs/match_barcodes/{chunk}.log'
    # Here we use process substitution like >() and <() to use zstandard to 
    # compress + decompress our inputs/outputs
    shell: "python3 {params.script} "
        " --R1_in <(zstd -dc {params.R1_in}) --R2_in <(zstd -dc {params.R2_in}) "
        " --R1_RNA >(zstd --fast=1 -q -o {output.R1_RNA}) "
        " --R2_RNA >(zstd --fast=1 -q -o {output.R2_RNA}) "
        " --R1_ATAC >(zstd --fast=1 -q -o {output.R1_ATAC}) "
        " --R2_ATAC >(zstd --fast=1 -q -o {output.R2_ATAC}) "
        " --ATAC_I2 {params.ATAC_I2} --RNA_I2 {params.RNA_I2} "
        " --BC1 {params.BC1} --BC2 {params.BC2} --BC3 {params.BC3} "
        " --json_stats {output.stats} "
        " 2> {log} "

# Remove adapter ends from the raw fastq reads
rule trim_adapters:
    input: 
        R1 = "{sample}/fastq/match_barcodes/{assay}/{chunk}_R1.fastq.zst",
        R2 = "{sample}/fastq/match_barcodes/{assay}/{chunk}_R2.fastq.zst",
    output:
        interleaved = temp("{sample}/fastq/trim_adapters/{assay}/{chunk}_interleaved.fastq.zst"),
        report_json = "{sample}/qc_stats/chunks/trim_adapters/{assay}/{chunk}.json"
    threads: 4
    log: '{sample}/logs/trim_adapters/{assay}/{chunk}.log'
    shell: "fastp --in1 {input.R1} --in2 {input.R2} "
        " -j {output.report_json} -G -Q -L -w {threads} 2> {log} "
        " --stdout | zstd --fast=1 -q -o {output.interleaved}"

# Alternative trim command
# "SeqPurge -a1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a2 CTGTCTCTTATACACATCTGACGCTGCCGACGA "
#     " -qcut 0 -ncut 0 "
#     " -threads {threads} -out1 {output.R1} -out2 {output.R2} "
#     " -in1 {input.R1} -in2 {input.R2} > {log}"


#fastp --adapter_sequence=GACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG
#TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
#AGCAGCCGTCGCAGTCTACACATATTCTCTGTC

#CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
#GACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG

# fastp alternative
# "fastp -i {input.fastq1_bc} -I {input.fastq2_bc} -o {output.fastq1_trim} -O {output.fastq2_trim}"
#        " -h {log.html} -j {log.json} -G -Q -L -w {threads} 2> {output.stats}"
