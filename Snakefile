configfile: "config.yaml"

rule all:
    input:
        config['output']

rule ga_align:
    input:
        gfa_input = config['graph_input'],
        reads = config['reads_file']
    output:
        config['output']
    params:
        cmd = "GraphAligner",
        x = "vg"
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x}"
