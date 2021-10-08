import os

def get_read_file_name():
    return os.path.basename(config['reads_file']).split('.')[0]

def get_read_file_path_wo_extension():
    index = config['reads_file'].rfind(".")
    return config['reads_file'][:index]

def get_fastq_for_sample(wildcards):
    return config[wildcards.sample]['fastq']

def get_config_for_sample(wildcards):
    return config[wildcards.sample]['config']

def get_ref(wildcards):
    return config['reference']

configfile: "config.yaml"

rule all:
    input:
        expand("{sample}.ga.tsv", sample=config['samples'])

rule gfa_gen:
    input:
        ref_file = get_ref,
        config_file = get_config_for_sample
    output:
        "{sample}.reference.gfa"
    params:
        cmd = "python",
        script = config['scripts_dir'] + "genome_str_graph_generator.py"
    conda: "ga.yaml"
    shell: 
        "{params.cmd} {params.script} --ref {input.ref_file} --config {input.config_file} > {output}"

rule ga_align:
    input:
        gfa_input = "{sample}.reference.gfa",
        reads = get_fastq_for_sample
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
