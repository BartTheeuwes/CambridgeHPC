#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH --time=1-11:59:59
#SBATCH -e cellranger.%N.%j.err
#SBATCH --output cellranger-log-%J.txt

export PATH=/rds/project/rds-SDzz0CATGms/users/bt392/software/cellranger-arc-2.0.1/cellranger-arc-2.0.1:$PATH

cellranger-arc count --id=$1 \
                     --reference=/rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                     --libraries=/rds/project/rds-SDzz0CATGms/users/bt392/10_Eomes_invitro_gut/original/$1.csv \
                     --localcores=32 \
                     --localmem=40