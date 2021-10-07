configfile: config.yaml

rule all:
    input:
        config[output]

rule ga_align:
    input:
        gfa_input = config[graph_input] 
        reads = config[reads]
    output:
        config[output]
    params:
        cmd="GraphAligner"
    conda: "str.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x vg"
