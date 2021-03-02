# Detect interesting SARS-CoV-2 spike protein mutations from Sanger sequencing data

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/covid-spike-classification/README.html)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/kblin/covid-spike-classification?style=flat)](https://hub.docker.com/r/kblin/covid-spike-classification)

`covid-spike-classification` is a script to call interesting SARS-CoV-2 spike protein mutations
from Sanger sequencing to support the Danish COVID-19 monitoring efforts.

Using Sanger-sequenced RT-PCR product of the spike protein, this tool should pick up all relevant
mutations currently tracked (see [`covid_spike_classification/core.p`](https://github.com/kblin/covid-spike-classification/blob/main/covid_spike_classification/core.py#L15-L35)
for the full list of tracked mutations) and give a table with one row per sample and a
yes/no/failed column per tracked mutation.

This workflow is built and maintained at https://github.com/kblin/covid-spike-classification

## Installation

`covid-spike-classification` is distributed via this git repository, pypi or bioconda.


### Bioconda

Installing via bioconda is the fastest way to get up and running:

```sh
conda create -n csc -c conda-forge -c bioconda covid-spike-classification
conda activate csc
```

### git & pypi


When installing via git or pypi, you first need to install the external binary dependencies.


`covid-spike-classification` depends on three excellent tools to do most of the work:

* tracy (versions 0.5.3 & 0.5.7 tested)
* bowtie2 (version 2.4.2 tested)
* samtools (versions 1.10 & 1.11 tested)

If you have `conda` installed, the easiest way to get started is to just install these via calling
```sh
git clone https://github.com/kblin/covid-spike-classification.git
cd covid-spike-classification
conda env create -n csc -f environment.yml
conda activate csc
pip install .
```

### Docker, Podman, Singularity

While not technically an installation method, `covid-spike-classification` is also shipped as an OCI container.
To use it, you ideally run the container from a workflow management system like [Snakemake](https://snakemake.github.io/)
or [Nextflow](https://www.nextflow.io/) that will take care of mounting filesystems into the container for you.

The OCI container image is available from the Docker Hub [`kblin/covid-spike-classification`](https://hub.docker.com/r/kblin/covid-spike-classification)
repository.


## Setup

You also need to generate the samtools and bowtie2 indices for your reference genome. We ship a
copy of NC\_045512 and a script to generate these indices:

```sh
conda activate csc
cd ref
./build_indices.sh
cd ..
```

## Usage

Assuming you used above instructions to install via conda, you can run the tool like this:

```sh
conda activate csc
covid-spike-classification --reference /path/to/your/reference.fasta --outdir /path/to/result/dir /path/to/sanger/reads/dir_or.zip
```

Notably, you can provide the input either as a ZIP file or as a directory, as long as they directly contain the ab1 files you want
to run the analysis on.

See also the `--help` output for more detailed usage information.


## License
All code is available under the Apache License version 2, see the
[`LICENSE`](LICENSE) file for details.
