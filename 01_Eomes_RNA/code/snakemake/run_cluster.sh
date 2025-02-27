#!/bin/bash
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time 4:00:00
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt
snakemake --cores -j 99 --latency-wait 90 -p --cluster "sbatch -p skylake-himem -A gottgens-sl2-cpu -n {threads} --time 11:00:00 --mem {resources.mem_mb}M"