#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -pe threads 4
#$ -V

snakemake --cores 4 --keep-going --latency-wait 30 > logs/$(date "+%Y_%m_%d %H:%M:%S").log

