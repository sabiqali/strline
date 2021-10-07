import os

def get_read_file_name():
    return os.path.basename(config['reads_file']).split('.')[0]

def get_read_file_path_wo_extension():
    index = config['reads_file'].rfind(".")
    return config['reads_file'][:index]

configfile: "config.yaml"

rule all:
    input:
        config['output_dir'] + get_read_file_name() + "_ga.tsv"

rule ga_align:
    input:
        gfa_input = config['graph_input'],
        reads = config['reads_file']
    output:
        config['output_dir'] + get_read_file_name() + ".gaf"
    params:
        cmd = "GraphAligner",
        x = "vg"
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x}"

rule ga_counter:
    input:
        gaf_input = config['output_dir'] + get_read_file_name() + ".gaf",
        name = get_read_file_name()
    output:
        config['output_dir'] + get_file_name(config['reads_file']) + "_ga.tsv"
    params:
        cmd = 'python',
        script = 'graph_counter.py'
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --name {input.gaf_input} --name {input.name} > {output}"