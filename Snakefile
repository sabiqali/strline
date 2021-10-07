import os

def get_file_name(file_path):
    return os.path.basename(file_path).split('.')[0]

def get_file_path_wo_extension(file_path):
    index = file_path.rfind(".")
    return file_path[:index]

configfile: "config.yaml"

rule all:
    input:
        config['output_dir'] + get_file_name(config['reads_file']) + "_ga.tsv"

rule ga_align:
    input:
        gfa_input = config['graph_input'],
        reads = config['reads_file']
    output:
        config['output_dir'] + get_file_name(config['reads_file']) + ".gaf"
    params:
        cmd = "GraphAligner",
        x = "vg"
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x}"

rule ga_counter:
    input:
        gaf_input = config['output_dir'] + get_file_name(config['reads_file']) + ".gaf",
        name = get_file_name(config['reads_file'])
    output:
        config['output_dir'] + get_file_name(config['reads_file']) + "_ga.tsv"
    params:
        cmd = 'python',
        script = 'graph_counter.py'
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --name {input.gaf_input} --name {input.name} > {output}"