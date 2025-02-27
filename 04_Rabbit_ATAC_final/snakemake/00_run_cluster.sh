#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time 23:59:59
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt
snakemake --cores 12 -j 99 --latency-wait 90 -p --cluster "sbatch -p skylake-himem -A gottgens-sl2-cpu -n {threads} --time 23:59:59 --mem {resources.mem_mb}M"