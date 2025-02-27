#!/bin/bash
#SBATCH -p sapphire
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time 4:00:00
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt
snakemake --cores 5 -j 10 --latency-wait 90 -p --cluster "sbatch -p sapphire -A gottgens-sl2-cpu -n {threads} --time 11:00:00 --mem {resources.mem_mb}M"
