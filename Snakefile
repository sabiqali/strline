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

def get_strscore_config_for_sample(wildcards):
    return config[wildcards.sample]['strscore_config']

def get_ref(wildcards):
    return config['reference']

def get_output_dir(wildcards):
    return config['output_dir']

def get_strique_index_for_sample(wildcards):
    return config[wildcards.sample]['strique_index']

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
        script = config['scripts_dir'] + "genome_str_graph_generator.py",
        memory_per_thread="1G"
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
        x = "vg",
        memory_per_thread="1G"
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x} --multiseed-DP 1"

rule ga_counter:
    input:
        gaf_input = "graphaligner/{sample}.gaf"
    output:
        "{sample}.ga.tsv"
    params:
        cmd = "python",
        script = config['scripts_dir'] + "parse_gaf.py",
        memory_per_thread="1G"
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --input {input.gaf_input} > {output}"

rule strscore_count_split:
    input:
        bam_file="splits/{sample}.split{splitID}.sorted.bam",
        bai_file="splits/{sample}.split{splitID}.sorted.bam.bai",
        reads_file="splits/{sample}.split{splitID}.fastq",
        ref_file = get_ref,
        config_file = get_strscore_config_for_sample
    output:
        "strscore/{sample}.split{splitID}_strscore.tsv"
    threads: 8
    params:
        cmd = "python",
        script = config['scripts_dir'] + "strscore_plasmids.py",
        memory_per_thread="2G"
    shell:
        "{params.cmd} {params.script} --bam {input.bam_file} --read {input.reads_file} --ref {input.ref_file} --config {input.config_file}> {output}"

rule strscore_merge:
    input:
        dynamic("strscore/{sample}.split{splitID}_strscore.tsv")
    output:
        "{sample}.strscore.tsv"
    threads: 1
    params:
        memory_per_thread="1G"
    shell:
        "cat {input} | awk 'NR == 1 || $0 !~ /strand/' > {output}"

rule compile_reads:
    input:
        ga_out = "{sample}.ga.tsv",
        strscore_out = "{sample}.strscore.tsv",
        strique_out = "{sample}.strique.tsv"
    output:
        "{sample}.compiled.tsv"
    threads: 1
    params:
        cmd = "python",
        script = config['scripts_dir'] + "merge_method_calls.py",
        memory_per_thread="1G"
    conda: "ga.yaml"
    shell:
        "{params.cmd} {params.script} --graphaligner {input.ga_out} --strscore {input.strscore_out} --strique {input.strique_out} > {output}"

rule split_index:
    input:
        strique_index = get_strique_index_for_sample
    output:
        dynamic("splits/{sample}.index.split{splitID}")
    threads: 1
    params:
        memory_per_thread="1G"
    shell:
        "cut -f2 {input.strique_index} | split -l {config[strique_reads_per_chunk]} - splits/{wildcards.sample}.index.split"

rule subset_reads:
    input:
        chunk="splits/{sample}.index.split{splitID}",
        reads_file=get_fastq_for_sample
    output:
        "splits/{sample}.split{splitID}.fastq"
    threads: 1
    params:
        memory_per_thread="1G"
    shell:
        "seqtk subseq {input.reads_file} {input.chunk} > {output}"

rule map:
    input:
        reads="{prefix}.fastq",
        ref=get_ref
    output:
        "{prefix}.sorted.bam"
    threads: 8
    params:
        memory_per_thread="4G"
    shell:
         "minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} | samtools sort -T {wildcards.prefix} > {output}"

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads: 1
    params:
        memory_per_thread="2G"
    shell:
         "samtools index {input}"

rule strique:
    input:
        bam="splits/{sample}.split{splitID}.sorted.bam",
        index=get_strique_index_for_sample,
        config=get_config_for_sample
    output:
        "strique/{sample}.split{splitID}_strique.tsv"
    threads: 8
    params:
        memory_per_thread="2G"
    shell:
        "samtools view -F 4 {input.bam} | python3 {config[strique]} count {input.index} {config[strique_pore_model]} {input.config} --out {output} --t {threads}"

rule strique_merge:
    input:
        dynamic("strique/{sample}.split{splitID}_strique.tsv")
    output:
        "{sample}.strique.tsv"
    threads: 1
    params:
        memory_per_thread="1G"
    shell:
        "cat {input} | awk 'NR == 1 || $0 !~ /target/' > {output}"
