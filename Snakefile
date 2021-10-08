import os

def get_read_file_name():
    return os.path.basename(config['reads_file']).split('.')[0]

def get_read_file_path_wo_extension():
    index = config['reads_file'].rfind(".")
    return config['reads_file'][:index]

configfile: "config.yaml"

rule all:
    input:
        expand("{sample}.ga.tsv", sample=config['samples'])

rule ga_align:
    input:
        gfa_input = "{sample}.reference.gfa",
        reads = "{sample}.fastq"
    output:
        "graphaligner/{sample}.gaf"
    params:
        cmd = "GraphAligner",
        x = "vg"
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x}"

rule ga_counter:
    input:
        gaf_input = "graphaligner/{sample}.gaf"
    output:
        "{sample}.ga.tsv"
    params:
        cmd = "python",
        script = config['scripts_dir'] + "graph_counter.py"
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --input {input.gaf_input} > {output}"
