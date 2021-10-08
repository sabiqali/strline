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

def get_bam_for_sample(wildcards):
    return config[wildcards.sample]['bam']

def get_ref(wildcards):
    return config['reference']

def get_output_dir(wildcards):
    return config['output_dir']

configfile: "config.yaml"

rule all:
    input:
        expand("{sample}.compiled.tsv", sample=config['samples'])

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

rule strscore_count:
    input:
        bam_file = get_bam_for_sample
        reads_file = get_fastq_for_sample
        ref_file = get_ref
        config_file = get_config_for_sample
    output:
        "{sample}.strscore.gaf"
    params:
        cmd = "python",
        script = config['scripts_dir'] + "strscore_plasmids.py"
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --bam {input.bam_file} --read {input.reads_file} --ref {input.ref_file} --config {input.config_file}> {output}"

rule compile_reads:
    input:
        ga_out = "{sample}.ga.tsv",
        strscore_out = "{sample}.strscore.gaf"
    output:
        "{sample}.compiled.tsv"
    params:
        cmd = "python",
        script = config['scripts_dir'] + "merge_method_calls.py"
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --graphaligner {input.ga_out} --strscore {input.strscore_out} > {output}"