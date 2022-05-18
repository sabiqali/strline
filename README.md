# Strline
Software pipeline to analyse STR loci from long read data. Strline can count the number of repeats in a repeat expansion, as well as compare it to other state-of-the-art repeat counting packages, to give you an analyses of per basecaller accuracy. Or, simply give you the count in a tabular format for further downstream analysis.

![Strline-flowchart](https://user-images.githubusercontent.com/39552869/168863661-9b8a9cd0-1a53-40a2-a805-82aeaed1a5bc.jpg)

## Release notes
* 0.1.0: Initial release with all the methods operating on all the all the different basecallers. Different methods and basecallers can be toggled on or off depending on need.

## Dependencies
For all packages apart from STRique Python 3.6 and above. Developed and tested on Python 3.7.10. Dependencies include:

* [STRique](https://github.com/giesselmann/STRique)
* [straglr](https://github.com/bcgsc/straglr)
* [tandem-genotyopes](https://github.com/mcfrith/tandem-genotypes)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://github.com/samtools/samtools)
* [seqtk](https://github.com/lh3/seqtk)
* [trf](https://github.com/Benson-Genomics-Lab/TRF)
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [anaconda/miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Installation instructions

### Installing Dependency: STRique

Please follow instructions [here](https://strique.readthedocs.io/en/latest/installation/prerequisites/) to install `STRique`. You will need to create a python virtual environment and have it activated when using `Strline`.

### Installing Dependency: straglr

Please make sure that you have `trf` and `blast` in your `$PATH`, to make sure that straglr works in the pipeline.

### Installing Strline

You can download latest code from github and install all dependencies except aforementioned ones as follows:

```
#clones the repository
git clone --recursive https://github.com/sabiqali/strline.git

#makes sure that the submodules have been initialised
git submodule update --init --recursive

cd strline

#install required packages and create a conda env
conda env create -n strline --file strline.yml
```

## Config file creation

The config file template should be downloaded along with the repository. Please open this file and make the required changes to the config file before running the pipeline. The config file has the instructions in it as to what to change. The pipeline will not run without these changes. 

## Running the pipeline

### Running the pipeline singularly

Copy the `config.yaml` file to the directory that you want to run the workflow. If you want to basecall the reads, please make sure you are on a computer with the GPU accessible and run the following:

```
snakemake -s /path/to/snakefile --rerun-incomplete --keep-going --latency-wait 60 --cores <specify_number_cores(1 if unsure)> plots
```

### Running the pipeline on a cluster/grid engine

Copy the `config.yaml` file to the directory that you want to run the workflow. There are a few different grid engines, so the exact format to run the workflow may be different for your particular grid engine:

```
snakemake --rerun-incomplete -s /path/to/snakefile --keep-going --jobs 500 --latency-wait 120 --cluster "qsub -cwd -V -o snakemake_all.output.log -e snakemake_all.error.log -N {rule} -pe smp {threads} -l h_vmem={params.memory_per_thread} {params.extra_cluster_opt} -l h_stack=32M -P <project_name> -b y" plots
```

You will have to replace `queue_name` and `project_name` with the necessary values to run on your cluster. `queue_name` is located inside the `Snakefile`.
