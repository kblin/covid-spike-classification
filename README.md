# Detect interesting SARS-CoV-2 spike protein variants from Sanger sequencing data

`covid-spike-classification` is a script to call interesting SARS-CoV-2 spike protein variants
from Sanger sequencing to support the Danish COVID-19 monitoring efforts.

Using Sanger-sequenced RT-PCR product of the spike protein, this tool should pick up all relevant
mutations currently tracked and give a table with one row per sample and a yes/no/failed column per
tracked mutation.

This workflow is built and maintained at https://github.com/kblin/covid-spike-classification

## Installation

For now, `covid-spike-classification` is distributed via this git repository.

`covid-spike-classification` depends on three excellent tools to do most of the work:

* tracy (version 0.5.3 tested)
* bowtie2 (version 2.4.2 tested)
* samtools (version 1.11 tested)

If you have `conda` installed, the easiest way to get started is to just install these via calling
```sh
git clone https://github.com/kblin/covid-spike-classification.git
cd covid-spike-classification
conda env create -n csc -f environment.yml
conda activate csc
pip install .
```

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
