#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -pe threads 10
#$ -V

snakemake --cores 10 --keep-going --latency-wait 30 > logs/$(date "+%Y_%m_%d %H:%M:%S").log

