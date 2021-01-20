#!/bin/bash

readonly REFERENCE="NC_045512.fasta"

samtools faidx ${REFERENCE}

bowtie2-build ${REFERENCE} ${REFERENCE%.*}.index
