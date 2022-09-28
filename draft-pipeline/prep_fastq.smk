
# sequencing:
#   220920:
#     type: bcl
#     run_dir: /oak/stanford/groups/wjg/seqruns/220920_A00509_0608_AHNT5JDSX3
#     ATAC_I2:
#       SL2_b1_lung: CTCTCTAT
#       SL4_b1_lung: AGAGTAGA
#     RNA_I2:
#       SL18_b1_lung: AGGTTGCC
#       SL20_b1_lung: CAGCAACG
#   20220726: 
#       type: fastq
#       ATAC_R1: /oak/stanford/groups/wjg/bliu/data/20220726_shareseq_rna_contam/fastq/raw/Undetermined_S0_R1_001.fastq.gz
#       ATAC_R2: /oak/stanford/groups/wjg/bliu/data/20220726_shareseq_rna_contam/fastq/raw/Undetermined_S0_R2_001.fastq.gz
#       RNA_R1: /oak/stanford/groups/wjg/bliu/data/20220726_shareseq_rna_contam/fastq/raw/Undetermined_S0_R1_001.fastq.gz
#       RNA_R2: /oak/stanford/groups/wjg/bliu/data/20220726_shareseq_rna_contam/fastq/raw/Undetermined_S0_R2_001.fastq.gz

workdir: config["output_dir"]

# Fix up any run_ids that are not strings for keys
for run_id, run in list(config["sequencing"].items()):
    if not isinstance(run_id, str):
        config["sequencing"][str(run_id)] = run
        del config["sequencing"][run_id]

# Need: list of sequencing input files 
def get_sequencing_paths():
    """Get a list of all paths for sublibraries for the current sample
    """
    sequencing_paths = []
    for run_id, run in config["sequencing"].items():
        assert run["type"] in ["bcl", "fastq"]
        if run["type"] == "bcl":
            sequencing_paths += [
                f"{run_id}/ATAC/{sublib_id}" for sublib_id in run["ATAC_I2"]
            ]
            sequencing_paths += [
                f"{run_id}/RNA/{sublib_id}" for sublib_id in run["RNA_I2"]
            ]
        elif run["type"] == "fastq":
            sequencing_paths += [
                f"{run_id}/ATAC", f"{run_id}/RNA"
            ]
        else:
            assert False
    return sequencing_paths

def get_fastq_path(sequencing_path, read):
    """Take a sublibrary path and return the path to its R1 or R2 fastq"""
    sequencing_path = sequencing_path.split("/") 
    run_id, assay_type = sequencing_path[:2]
    if config["sequencing"][run_id]["type"] == "bcl":
        sublib_id = sequencing_path[2]
        return f"bcl2fastq/{run_id}/{assay_type}/{sublib_id}_{read}.fastq.gz"
    elif config["sequencing"][run_id]["type"] == "fastq":
        return config["sequencing"][run_id][assay_type.upper() + "_" + read]
    else:
        assert False

rule all:
    input: 
        expand("fastq/{sequencing_path}/read_count.txt", sequencing_path=get_sequencing_paths())

rule bcl2fastq:
    input:
        run_dir = lambda w: config["sequencing"][w.sequencing_run]["run_dir"]
    output:
        samplesheet = "bcl2fastq/{sequencing_run}_SampleSheet.csv",
        results_dir = temp(directory("bcl2fastq/raw/{sequencing_run}/")),
    params:
        cluster_memory = "32G",
        cluster_time = "10:00:00"
    threads: 20
    log: "logs/bcl2fastq/{sequencing_run}.log"
    benchmark: "logs/bcl2fastq/{sequencing_run}.runtime"
    script: "scripts/run_bcl2fastq.py"


def get_fastq_input(sequencing_path):
    """Take a sublibrary path and return a list of input dependencies"""
    sequencing_path = sequencing_path.split("/") 
    run_id, assay_type = sequencing_path[:2]
    if config["sequencing"][run_id]["type"] == "bcl":
        sublib_id = sequencing_path[2]
        return f"bcl2fastq/raw/{run_id}/"
    elif config["sequencing"][run_id]["type"] == "fastq":
        return [
            config["sequencing"][run_id][assay_type.upper() + "_R1"],
            config["sequencing"][run_id][assay_type.upper() + "_R2"],
        ]
    else:
        assert False

rule bcl2fastq_index_to_read_names:
    input: lambda w: get_fastq_input(w.sequencing_path)
    output:
        R1 = "bcl2fastq/{sequencing_path}_R1.fastq.gz",
        R2 = "bcl2fastq/{sequencing_path}_R2.fastq.gz",
    params: 
        R1 = "bcl2fastq/raw/{sequencing_path}_R1.fastq.gz",
        R2 = "bcl2fastq/raw/{sequencing_path}_R2.fastq.gz",
        I1 = "bcl2fastq/raw/{sequencing_path}_I1.fastq.gz",
        I2 = "bcl2fastq/raw/{sequencing_path}_I2.fastq.gz",
        script = srcdir("scripts/index_fastqs_to_readnames.py")
    threads: 3
    shell: "python {params.script} "
            " --R1 {params.R1} --R2 {params.R2} --I1 {params.I1} --I2 {params.I2} "
            " --R1_out {output.R1} --R2_out {output.R2}"
    
rule count_reads:
    input: lambda w: get_fastq_path(w.sequencing_path, "R1")
    output: "fastq/{sequencing_path}/read_count.txt"
    shell: "gunzip -c {input} | awk -c 'END{{print int(NR/4)}}' > {output}"


