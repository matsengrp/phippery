#!/bin/bash

set -e
#source /app/lmod/lmod/init/profile
##
#module load nextflow
#module load Singularity
#export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C phipflow.config.docker \
    run phip-flow/PhIP-Flow.nf
