#! /bin/bash
snakemake --log-handler-script log.py --resources mem_mb=300000 -j 90
