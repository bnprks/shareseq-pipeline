# Fastq preparation to run before main shareseq pipeline
# Author: Ben Parks
# Last Modified: 9/28/22

# 1. Optionally demultiplex a run using bcl2fastq
# 2. Count the number of reads per fastq file

workdir: config["output_dir"]

wildcard_constraints:
    sequencing_path = "(ATAC|RNA)/([^/]+/)?[^/]+", # Sequencing path is 2-3 folders

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

localrules: all
rule all:
    input: 
        expand("{sequencing_path}/read_count.txt", sequencing_path=get_sequencing_paths("ATAC") + get_sequencing_paths("RNA"))

#############################
### Direct BCL extraction
#############################

localrules: bcl2fastq_samples
rule bcl2fastq_samples:
    output: 
        samples = "bcl2fastq/samples_{sequencing_run}.csv"
    script: "scripts/prep_fastq/prep_samplesheet.py"

rule bcl2fastq:
    input:
        run_dir = lambda w: config["sequencing"][w.sequencing_run]["run_dir"],
        samples = rules.bcl2fastq_samples.output,
    output:
        results_dir = temp(directory("bcl2fastq/raw/{sequencing_run}/")),
    params:
        script = srcdir("scripts/prep_fastq/run_bcl2fastq.py"),
    resources:
        runtime = 8 * 60,
        mem_mb = 32_000,
    threads: 24
    log: "bcl2fastq/logs/{sequencing_run}.log"
    benchmark: "logs/bcl2fastq/{sequencing_run}.runtime"
    shell: "python {params.script}"
        " --samples {input.samples} "
        " --input {input.run_dir} "
        " --output {output.results_dir} "
        " --threads {threads} "
        " --log {log}"

def bcl2fastq_dependency(sequencing_path):
    """Take a sublibrary path and return a list of input dependencies"""
    run_id = sequencing_path.split("/")[1]
    assert config["sequencing"][run_id]["type"] == "bcl"
    return f"bcl2fastq/raw/{run_id}/"

def fastq_path(sequencing_path, read):
    """Take a sublibrary path and return the path to its R1 or R2 fastq"""
    assay_type, run_id = sequencing_path.split("/")[:2] 
    if config["sequencing"][run_id]["type"] == "bcl":
        return f"bcl2fastq/{sequencing_path}_{read}.fastq.gz"
    elif config["sequencing"][run_id]["type"] == "fastq":
        return config["sequencing"][run_id][assay_type.upper() + "_" + read]
    else:
        assert False

def bcl2fastq_output(sequencing_path, read):
    """Take a bcl sublibrary path and return the path to its bcl2fastq path"""
    assay_type, run_id = sequencing_path.split("/")[:2]
    assert config["sequencing"][run_id]["type"] == "bcl"
    sublib_id = sequencing_path.split("/")[2]
    return f"bcl2fastq/raw/{run_id}/{assay_type}_{sublib_id}_{read}.fastq.gz"

rule bcl2fastq_index_to_read_names:
    input: lambda w: bcl2fastq_dependency(w.sequencing_path)
    output:
        R1 = "bcl2fastq/{sequencing_path}_R1.fastq.gz",
        R2 = "bcl2fastq/{sequencing_path}_R2.fastq.gz",
    params: 
        R1 = lambda w: bcl2fastq_output(w.sequencing_path, "R1"),
        R2 = lambda w: bcl2fastq_output(w.sequencing_path, "R2"),
        I1 = lambda w: bcl2fastq_output(w.sequencing_path, "I1"),
        I2 = lambda w: bcl2fastq_output(w.sequencing_path, "I2"),
        script = srcdir("scripts/prep_fastq/fastq_index_to_readname.py")
    threads: 3
    shell: "python {params.script} "
            " --R1 {params.R1} --R2 {params.R2} --I1 {params.I1} --I2 {params.I2} "
            " --R1_out {output.R1} --R2_out {output.R2}"
    
#############################
### Count reads
#############################

rule count_reads:
    input: lambda w: fastq_path(w.sequencing_path, "R1")
    output: "{sequencing_path}/read_count.txt"
    shell: "gunzip -c {input} | awk -c 'END{{print int(NR/4)}}' > {output}"


