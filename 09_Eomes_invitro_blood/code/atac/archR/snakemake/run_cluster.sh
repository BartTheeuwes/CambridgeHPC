#!/bin/bash
#SBATCH -p cclake # skylake-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time 11:00:00
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt
snakemake --cores 16 -j 99 -r --latency-wait 90 -p --cluster "sbatch -p cclake -A gottgens-sl2-cpu -n {threads} --time 11:00:00 --mem {resources.mem_mb}M" #skylake-himem