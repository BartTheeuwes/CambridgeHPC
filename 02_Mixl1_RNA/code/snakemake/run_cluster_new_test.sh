#!/bin/bash
#SBATCH -p cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time 12:00:00
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-ccl 

snakemake --cores -j 99 --latency-wait 90 -p --cluster "sbatch -p cclake-himem -A gottgens-sl2-cpu -n {threads} --time 11:00:00 --mem {resources.mem_mb}M"

# Like this works from command line it seems: 
#  snakemake --use-conda --cores 4 -j 2 --latency-wait 90 -p --cluster "sbatch -p sapphire -A gottgens-sl2-cpu -n {threads} --time 11:00:00 --mem {resources.mem_mb}M"