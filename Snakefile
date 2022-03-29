import os

def get_fastq_for_sample(wildcards):
    return "fastq/" + wildcards.sample + "." + wildcards.basecall_config + ".fastq"

def get_raw_file_input_path(wildcards):
    return config['fast5']

def get_fast5_path(wildcards):
    return config[wildcards.data_type]['fast5']

def get_fast5_path_for_sample(wildcards):
    dt = config[wildcards.sample]['data_type']
    return config[dt]['fast5']

def get_guppy_basecaller(wildcards):
    return config[wildcards.basecall_config]['path'] + "/guppy_basecaller"

def get_guppy_barcoder(wildcards):
    return config[wildcards.basecall_config]['path'] + "/guppy_barcoder"

def get_barcoding_kit(wildcards):
    return config[wildcards.data_type]['barcoding_kit']

def get_guppy_mode(wildcards):
    return config[wildcards.basecall_config]['mode']

def get_guppy_extra_opt(wildcards):
    return config[wildcards.basecall_config]['guppy_extra']

def get_basecalled_dir(wildcards):
    dt = config[wildcards.sample]['data_type']
    return "basecalled." + wildcards.basecall_config + "." + dt + "/"

def get_barcoded_dir(wildcards):
    dt = config[wildcards.sample]['data_type']
    return "barcoded." + wildcards.basecall_config + "." + dt + "/"

def get_bonito_exec(wildcards):
    return config['bonito']

#
# to support both singleplex and multiplex experiments this needs to be split
# into a function that gets the root dir, and one that gets the subdir to make
# snakemake happy and only run demultiplexing once.
#
def get_basecalled_root_dir_for_sample(wildcards):
    dt = config[wildcards.sample]['data_type']
    bk = config[dt]['barcoding_kit']
    bid = config[wildcards.sample]['barcode']
    # if barcoding is not enabled, merged all fastqs in the basecalled directory
    if bk == "none":
        p = get_basecalled_dir(wildcards) + "basecalled.fastq"
    else:
        # barcoding enabled
        p = get_barcoded_dir(wildcards) 
        #p = p + "barcode" + bid + ".fastq"
    return p

def get_basecalled_subdir_for_sample(wildcards):
    bc = get_barcode_id_for_sample(wildcards)
    if bc == "none":
        return ""
    else:
        return "barcode" + bc

def get_filename_if_barcoded(wildcards):
    dt = config[wildcards.sample]['data_type']
    bk = config[dt]['barcoding_kit']
    bid = config[wildcards.sample]['barcode']
    if bk == "none":
        p = ""
    else:
        # barcoding enabled
        #p = get_barcoded_dir(wildcards)
        p = "/barcode" + bid + ".fastq"
    return p

def get_sequencing_summary_for_sample(wildcards):
    return get_basecalled_dir(wildcards) + "basecalled_summary.tsv"

def get_barcode_id_for_sample(wildcards):
    return config[wildcards.sample]['barcode']

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

def get_microsat_file(wildcards):
    #return config['microsat']
    return get_data_type_config_for_sample(wildcards.sample)['microsat']

def get_copy_number(wildcards):
    #return config['cn']
    return get_data_type_config_for_sample(wildcards.sample)['cn']

def get_graphaligner_mode_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['graphaligner_mode']

def get_methods_for_sample(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['methods']

def get_compile_input_for_sample(wildcards):
    methods = get_data_type_config_for_sample(wildcards.sample)['methods']
    return expand("{sample}.{basecall_config}.{method}.tsv", sample=wildcards.sample, basecall_config=wildcards.basecall_config, method=methods)

def get_maximum_length_for_plot(wildcards):
    return get_data_type_config_for_sample(wildcards.sample)['maximum_length_for_plot']

configfile: "config.yaml"

rule all:
    input:
        expand("{sample}.{basecall_config}.compiled_singletons_only.tsv", sample=config['samples'], basecall_config=config['basecall_configs'])

rule plots:
    input:
        expand("{sample}.compiled_singletons_only_distributions_by_basecaller.pdf", sample=config['samples'])

#
# Helpers
#
rule last_index:
    input:
        ref_file=get_ref_for_sample
    output:
        "{sample}_refdb.bck",
        "{sample}_refdb.des",
        "{sample}_refdb.prj",
        "{sample}_refdb.sds",
        "{sample}_refdb.ssp",
        "{sample}_refdb.suf",
        "{sample}_refdb.tis"
    threads: 1
    params:
        memory_per_thread="256G",
        extra_cluster_opt=""
    shell:
        "lastdb -P8 -uNEAR {wildcards.sample}_refdb {input.ref_file}"

rule last_sg_rates:
    input:
        "{sample}_refdb.bck",
        "{sample}_refdb.des",
        "{sample}_refdb.prj",
        "{sample}_refdb.sds",
        "{sample}_refdb.ssp",
        "{sample}_refdb.suf",
        "{sample}_refdb.tis",
        reads_file=get_fastq_for_sample
    output:
        "{sample}.{basecall_config}.par"
    threads: 1
    params:
        memory_per_thread="256G",
        extra_cluster_opt=""
    shell:
        "last-train -P8 -Q0 {wildcards.sample}_refdb {input.reads_file} > {output}"

rule last_align:
    input:
        "{sample}_refdb.bck",
        "{sample}_refdb.des",
        "{sample}_refdb.prj",
        "{sample}_refdb.sds",
        "{sample}_refdb.ssp",
        "{sample}_refdb.suf",
        "{sample}_refdb.tis",
        par_file="{sample}.{basecall_config}.par",
        reads_file=get_fastq_for_sample
    output:
        "{sample}.{basecall_config}.maf"
    threads: 1
    params:
        memory_per_thread="256G",
        extra_cluster_opt=""
    shell:
        "lastal -P8 -p {input.par_file} {wildcards.sample}_refdb {input.reads_file} | last-split > {output}"

rule split_index:
    input:
        index = get_sequencing_summary_for_sample
    output:
        dynamic("splits/{sample}.{basecall_config}.index.split{splitID}")
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        "cut -f2 {input.index} | grep -v read_id | split -l {config[strique_reads_per_chunk]} - splits/{wildcards.sample}.{wildcards.basecall_config}.index.split"

rule subset_reads:
    input:
        chunk="splits/{sample}.{basecall_config}.index.split{splitID}",
        reads_file=get_fastq_for_sample
    output:
        "splits/{sample}.{basecall_config}.split{splitID}.fastq"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        "seqtk subseq {input.reads_file} {input.chunk} > {output}"

rule map_split:
    input:
        reads="splits/{sample}.{basecall_config}.split{splitID}.fastq",
        ref=get_ref_for_sample
    output:
        "splits/{sample}.{basecall_config}.split{splitID}.sorted.bam",
    threads: 8
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
         "minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} | samtools sort -T {wildcards.sample} > {output}"

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads: 1
    params:
        memory_per_thread="2G",
        extra_cluster_opt=""
    shell:
         "samtools index {input}"
#
# Basecaller
#
#rule basecall_reads:
#	input:
#		fast5 = get_fast5_path
#	output:
#		dir=directory("basecalled.{basecall_config}.{data_type}"), ss="basecalled.{basecall_config}.{data_type}/sequencing_summary.txt"
#	threads: 8
#	params:
#		mode=get_guppy_mode,
#		guppy_location=get_guppy_basecaller,
#		guppy_extra_opt=get_guppy_extra_opt,
#		memory_per_thread="8G",
#		extra_cluster_opt="-q gpu.q -l gpu=2"
#	shell:
#		"{params.guppy_location} --num_callers 8 --input_path {input.fast5} --save_path {output.dir} -c dna_r9.4.1_450bps_{params.mode}.cfg {params.guppy_extra_opt} -x 'cuda:0 cuda:1'" 
#
#rule demux_reads:
#	input:
#		"basecalled.{basecall_config}.{data_type}"
#	output:
#		directory("barcoded.{basecall_config}.{data_type}")
#	threads: 8
#	params:
#		guppy_location=get_guppy_barcoder,
#		barcoding_kit=get_barcoding_kit,
#		memory_per_thread="1G",
#		extra_cluster_opt="-q gpu.q -l gpu=2"
#	shell:
#		"{params.guppy_location} --barcode_kits {params.barcoding_kit} --recursive -i {input} -s barcoded.{wildcards.basecall_config}.{wildcards.data_type} -x 'cuda:0 cuda:1'"
#
#rule merge_reads:
#    input:
#        dir=get_basecalled_root_dir_for_sample
#    output:
#        "fastq/{sample}.{basecall_config}.fastq"
#    threads: 1
#    params:
#        memory_per_thread="1G",
#        extra_cluster_opt="",
#        subdir = get_basecalled_subdir_for_sample
#    shell:
#        "find {input.dir}/{params.subdir} -name \"*.fastq\" -exec cat {{}} + > {output}"

#
# Bonito
#
rule bonito_basecall:
    input:
        fast5 = get_fast5_path
    output:
        fastq_out="basecalled.{basecall_config}.{data_type}/basecalled.fastq", ss="basecalled.{basecall_config}.{data_type}/basecalled_summary.tsv"
    threads: 8
    params:
        mode=get_guppy_mode,
        bonito_exec=get_bonito_exec,
        memory_per_thread="8G",
        extra_cluster_opt="-q gpu.q -l gpu=2"
    shell:
        "{params.bonito_exec} basecaller {params.mode} {input.fast5} > {output.fastq_out}"

rule bonito_demux:
    input:
        "basecalled.{basecall_config}.{data_type}/basecalled.fastq"
    output:
        directory("barcoded.{basecall_config}.{data_type}")
    threads: 8
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "qcat -f {input} -b {output}"

rule move_reads:
    input:
        get_basecalled_root_dir_for_sample
    output:
        "fastq/{sample}.{basecall_config}.fastq"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt="",
	bid=get_barcode_id_for_sample,
        bfile=get_filename_if_barcoded
    shell:
        "mv {input}{params.bfile} {output}" 
#
# Tandem Genotypes
#
rule tg_detect:
    input:
        alignment_file="{sample}.{basecall_config}.maf",
        microsat_file=get_microsat_file
    output:
        "{sample}.{basecall_config}.tg.tsv"
    threads:1
    params:
        copy_number=get_copy_number,
        script = srcdir("scripts/tandem-genotypes.py"),
        memory_per_thread="32G",
        extra_cluster_opt=""
    shell:
        "{params.script} {input.microsat_file} {input.alignment_file} --cn {params.copy_number} > {output}"

#rule tg_parse:
#    input:
#        tg_out="tandem_genotypes/{sample}.{basecall_config}.txt",
#        config=get_repeat_config_for_sample,
#        copy_number=get_copy_number_for_sample
#    output:
#        "{sample}.{basecall_config}.tg.tsv"
#    threads:1
#    params:
#        script = srcdir("scripts/parse_tandem_genotypes.py"),
#        memory_per_thread="1G",
#        extra_cluster_opt=""
#    shell:
#        "{params.script} --input {input.tg_out} --config {input.config} --copy_number {input.copy_number} > {output}"

#
# GraphAligner
#
rule gfa_gen:
    input:
        ref_file = get_ref_for_sample,
        config_file = get_repeat_config_for_sample
    output:
        "{sample}.reference.gfa"
    threads: 1
    params:
        script = srcdir("scripts/genome_str_graph_generator.py"),
        memory_per_thread="32G",
        extra_cluster_opt=""
    conda: "ga.yaml"
    shell: 
        "{params.script} --ref {input.ref_file} --config {input.config_file} > {output}"

rule ga_align:
    input:
        gfa_input = "{sample}.reference.gfa",
        reads = get_fastq_for_sample
    output:
        "graphaligner/{sample}.{basecall_config}.gaf"
    threads: 1
    params:
        cmd = "GraphAligner",
        x = "vg",
        memory_per_thread="256G",
        extra_cluster_opt="",
        graphaligner_mode=get_graphaligner_mode_for_sample
    conda: "ga.yaml"
    shell:
        "{params.cmd} -g {input.gfa_input} -f {input.reads} -a {output} -x {params.x} {params.graphaligner_mode}"

rule ga_counter:
    input:
        gaf_input = "graphaligner/{sample}.{basecall_config}.gaf"
    output:
        "{sample}.{basecall_config}.ga.tsv"
    params:
        cmd = "python",
        script = srcdir("scripts/parse_gaf.py"),
        memory_per_thread="1G",
        extra_cluster_opt=""
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
        memory_per_thread="2G",
        extra_cluster_opt=""
    shell:
        "{params.script} --bam {input.bam_file} --read {input.reads_file} --ref {input.ref_file} --config {input.config_file}> {output}"

rule strscore_merge:
    input:
        dynamic("strscore/{sample}.split{splitID}_strscore.tsv")
    output:
        "{sample}.strscore.tsv"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        "cat {input} | awk 'NR == 1 || $0 !~ /strand/' > {output}"

#
# Strique
#
rule strique_index:
    input:
        ss=get_sequencing_summary_for_sample,
        fast5_dir=get_fast5_path_for_sample
    output:
        "strique/{sample}.{basecall_config}.index.fofn"
    params:
        memory_per_thread="2G",
        extra_cluster_opt="",
        script=srcdir("scripts/generate_strique_index.py")
    shell:
        "python {params.script} --fast5_dir {input.fast5_dir} --summary {input.ss} > {output}"

rule strique:
    input:
        bam="splits/{sample}.{basecall_config}.split{splitID}.sorted.bam",
        index="strique/{sample}.{basecall_config}.index.fofn",
        config=get_repeat_config_for_sample
    output:
        "strique/{sample}.{basecall_config}.split{splitID}_strique.tsv"
    threads: 8
    params:
        memory_per_thread="2G",
        extra_cluster_opt=""
    shell:
        "samtools view -F 4 {input.bam} | python3 {config[strique]} count {input.index} {config[strique_pore_model]} {input.config} --out {output} --t {threads}"

rule strique_merge:
    input:
        dynamic("strique/{sample}.{basecall_config}.split{splitID}_strique.tsv")
    output:
        "{sample}.{basecall_config}.strique.tsv"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        "cat {input} | awk 'NR == 1 || $0 !~ /target/' > {output}"

#
# simplecount
#
rule simplecount:
    input:
        reads=get_fastq_for_sample,
        config=get_repeat_config_for_sample
    output:
        "{sample}.{basecall_config}.simplecount.tsv"
    params:
        script = srcdir("scripts/simple_count.py"),
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        "{params.script} --config {input.config} {input.reads} > {output}"

#
# Merge/analysis
#
rule compile_results:
    input:
        get_compile_input_for_sample
    output:
        "{sample}.{basecall_config}.compiled.tsv"
    threads: 1
    params:
        script = srcdir("scripts/merge_method_calls.py"),
        memory_per_thread="1G",
        extra_cluster_opt=""
    conda: "ga.yaml"
    shell:
        "{params.script} {input} > {output}"

rule get_singletons:
    input:
        "{sample}.{basecall_config}.ga.tsv"
    output:
        "{sample}.{basecall_config}.ga_singletons.txt"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        """grep -v read_name {input} | cut -f1 | sort | uniq -c | awk "{{ if(\$1 == 1) {{ print \$2 }} }}" > {output}"""

rule filter_results_singletons:
    input:
        all_results = "{sample}.{basecall_config}.compiled.tsv",
        singleton_ids = "{sample}.{basecall_config}.ga_singletons.txt"
    output:
        "{sample}.{basecall_config}.compiled_singletons_only.tsv"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt=""
    shell:
        """
        head -1 {input.all_results} > {output}
        grep -f {input.singleton_ids} {input.all_results} >> {output}
        """

rule plot_distributions:
    input:
        "{sample}.{basecall_config}.compiled_singletons_only.tsv"
    output:
        "{sample}.{basecall_config}.compiled_singletons_only_distributions.pdf"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt="",
        script = srcdir("scripts/plot_str_distributions.R"),
        max_length=get_maximum_length_for_plot
    shell:
        "Rscript {params.script} --input {input} --output {output} --maximum-length {params.max_length}"

rule plot_distributions_by_basecaller:
    input:
        expand("{{sample}}.{basecall_config}.{{results_type}}.tsv", basecall_config=config['basecall_configs'])
    output:
        "{sample}.{results_type}_distributions_by_basecaller.pdf"
    threads: 1
    params:
        memory_per_thread="1G",
        extra_cluster_opt="",
        script = srcdir("scripts/plot_str_distributions_by_basecaller.R"),
        max_length=get_maximum_length_for_plot
    shell:
        "Rscript {params.script} --output {output} --maximum-length {params.max_length} {input}"
