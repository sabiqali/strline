import os

def get_fastq_for_sample(wildcards):
    return config[wildcards.sample]['fastq']

def get_data_type_config_for_sample(sample):
    return config[config[sample]['data_type']]

def get_ref_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['reference']

def get_repeat_config_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['repeat_config']

def get_strscore_config_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['strscore_config']

def get_strique_index_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['strique_index']

def get_graphaligner_mode_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['graphaligner_mode']

def get_methods_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['methods']

def get_compile_input_for_sample(wildcards):
    methods = get_data_type_config_for_sample(wildcards.sample)['methods']
    return expand("{sample}.{method}.tsv", sample=wildcards.sample, method=methods)

def get_maximum_length_for_plot(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['maximum_length_for_plot']

configfile: "config.yaml"

rule all:
    input:
        expand("{sample}.compiled_singletons_only.tsv", sample=config['samples'])

rule plots:
    input:
        expand("{sample}.compiled_singletons_only_distributions.pdf", sample=config['samples'])

#
# Helpers
#
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

rule map_split:
    input:
        reads="splits/{sample}.split{splitID}.fastq",
        ref=get_ref_for_sample
    output:
        "splits/{sample}.split{splitID}.sorted.bam",
    threads: 8
    params:
        memory_per_thread="4G"
    shell:
         "minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} | samtools sort -T {wildcards.sample} > {output}"

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

#
# GraphAligner
#
rule gfa_gen:
    input:
        ref_file = get_ref_for_sample,
        config_file = get_repeat_config_for_sample
    output:
        "{sample}.reference.gfa"
    params:
        script = srcdir("scripts/genome_str_graph_generator.py"),
        memory_per_thread="1G"
    conda: "ga.yaml"
    shell: 
        "{params.script} --ref {input.ref_file} --config {input.config_file} > {output}"

rule ga_align:
    input:
        gfa_input = "{sample}.reference.gfa",
        reads = get_fastq_for_sample
    output:
        "graphaligner/{sample}.gaf"
    params:
        cmd = "GraphAligner",
        x = "vg",
        memory_per_thread="1G",
        graphaligner_mode=get_graphaligner_mode_for_sample
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x} {params.graphaligner_mode}"

rule ga_counter:
    input:
        gaf_input = "graphaligner/{sample}.gaf"
    output:
        "{sample}.ga.tsv"
    params:
        cmd = "python",
        script = srcdir("scripts/parse_gaf.py"),
        memory_per_thread="1G"
    conda: "ga.yaml"
    shell:
        "{params.script} --input {input.gaf_input} > {output}"

#
# strscore
#
rule strscore_count_split:
    input:
        bam_file="splits/{sample}.split{splitID}.sorted.bam",
        bai_file="splits/{sample}.split{splitID}.sorted.bam.bai",
        reads_file="splits/{sample}.split{splitID}.fastq",
        ref_file = get_ref_for_sample,
        config_file = get_strscore_config_for_sample
    output:
        "strscore/{sample}.split{splitID}_strscore.tsv"
    threads: 8
    params:
        script = srcdir("scripts/strscore_plasmids.py"),
        memory_per_thread="2G"
    shell:
        "{params.script} --bam {input.bam_file} --read {input.reads_file} --ref {input.ref_file} --config {input.config_file}> {output}"

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

#
# Strique
#
rule strique:
    input:
        bam="splits/{sample}.split{splitID}.sorted.bam",
        index=get_strique_index_for_sample,
        config=get_repeat_config_for_sample
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

#
# Merge/analysis
#
rule compile_results:
    input:
        get_compile_input_for_sample
    output:
        "{sample}.compiled.tsv"
    threads: 1
    params:
        script = srcdir("scripts/merge_method_calls.py"),
        memory_per_thread="1G",
    conda: "ga.yaml"
    shell:
        "{params.script} {input} > {output}"

rule get_singletons:
    input:
        "{sample}.ga.tsv"
    output:
        "{sample}.ga_singletons.txt"
    threads: 1
    params:
        memory_per_thread="1G"
    shell:
        """grep -v read_name {input} | cut -f1 | sort | uniq -c | awk "{{ if(\$1 == 1) {{ print \$2 }} }}" > {output}"""

rule filter_results_singletons:
    input:
        all_results = "{sample}.compiled.tsv",
        singleton_ids = "{sample}.ga_singletons.txt"
    output:
        "{sample}.compiled_singletons_only.tsv"
    threads: 1
    params:
        memory_per_thread="1G"
    shell:
        """
        head -1 {input.all_results} > {output}
        grep -f {input.singleton_ids} {input.all_results} >> {output}
        """

rule plot_distributions:
    input:
        "{sample}.compiled_singletons_only.tsv"
    output:
        "{sample}.compiled_singletons_only_distributions.pdf"
    threads: 1
    params:
        memory_per_thread="1G",
        script = srcdir("scripts/plot_str_distributions.R"),
        max_length=get_maximum_length_for_plot
    shell:
        "Rscript {params.script} --input {input} --output {output} --maximum-length {params.max_length}"
