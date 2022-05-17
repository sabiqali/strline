# Strline
Software pipeline to analyse STR loci from long read data. Strline can count the number of repeats in a repeat expansion, as well as compare it to other state-of-the-art repeat counting packages, to give you an analyses of per basecaller accuracy. Or, simply give you the count in a tabular format for further downstream analysis.

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

Please follow instructions [here](https://strique.readthedocs.io/en/latest/installation/prerequisites/) to install STRique. You will need to create a python virtual environment and have it activated when using Strline.

### Installing Dependency: straglr

Please make sure that you have `trf` and `blast` in your `$PATH`, to make sure that straglr works in the pipeline.

### Installing Strline

You can download latest code from github and install all dependencies except aforementioned ones as follows:

```
git clone --recursive https://github.com/sabiqali/strline.git
cd strline
conda env create -n strline --file strline.yml
```

## Config file creation


