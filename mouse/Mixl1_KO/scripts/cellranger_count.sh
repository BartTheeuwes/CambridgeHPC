#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -p clincloud-himem
#SBATCH -A GOTTGENS-CCLD-SL2-CPU
#SBATCH --time=0-11:59:59
#SBATCH -e slurm.%N.%j.err

# -----------------------------------------------------------
#title           :cellranger_count.sh
#author          :hpb29
#date            :20180404
#version         :0.7    

#description     :This script submits 10x libraries via
#                 SLURM to 'cellranger count' and expects 
#                 to be fed four positional arguments:
#                 $1 - library/sample id
#                 $2 - path to 10x reference
#                 $3 - path to fastq files dir
#                 $4 - expected number of cells

#usage            :On this file you are only supposed to change
#                 the SLURM parameters (if necessary). For the
#                 cellranger parameters it is best to submit
#                 through the 'run_cellranger_count.sh' file.
#                 
# ------------------------------------------------------------

module load cellranger/2.1.1

echo "Started processing library $1 on `date`"

cellranger count --id=$1 \
                 --sample=$1 \
                 --transcriptome=$2 \
                 --fastqs=$3 
                 
echo "Finishing processing library $1 on `date`"
