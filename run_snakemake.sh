#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake -k -c 72 --use-singularity --singularity-args "-B /projects,/gsc,/home"  --rerun-incomplete
