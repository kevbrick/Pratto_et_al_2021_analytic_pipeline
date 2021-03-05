#!/bin/bash

nextflow run -c `pwd`/accessoryFiles/config/nextflow.config.nf \
             -profile singularity \
             `pwd`/Pratto_et_al_2021_analyticPipeline.nf \
             --outdir output \
             --projectdir `pwd`
