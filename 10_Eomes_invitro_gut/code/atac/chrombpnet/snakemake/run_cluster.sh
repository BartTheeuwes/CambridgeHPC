#!/bin/bash
#SBATCH -p sapphire # skylake-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --time 23:55:00
#SBATCH --job-name snakemake
#SBATCH --output logs/snakemake-log-%J.txt

#snakemake --cores 5 -j 99 -r --latency-wait 90 --use-conda -p --cluster "sbatch --partition={resources.partition} --account={resources.account} --time={resources.time} --nodes={resources.N} --ntasks={resources.n} {resources.slurm_extra}"

snakemake --cores 5 -j 99 -r --latency-wait 90 --keep-incomplete --ignore-incomplete --use-conda --conda-base-path "/home/bt392/miniconda3/" -p --cluster "sbatch --partition={resources.partition} --account={resources.account} --time={resources.time} --nodes={resources.N} --ntasks={resources.n} {resources.slurm_extra}"

# 

#snakemake --cores 5 -j 10 -r --latency-wait 90 --use-conda -p --slurm
