#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -p skylake-himem
#SBATCH -A gottgens-sl2-cpu
#SBATCH --time=0-20:59:59
#SBATCH -e cellranger.%N.%j.err
#SBATCH --output cellranger-log-%J.txt
                
                cellranger-arc count --id=SITTA9 \
                       --reference=/rds/project/rds-SDzz0CATGms/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                       --libraries=/rds/project/rds-SDzz0CATGms/users/bt392/10_Eomes_invitro_gut/data/raw/RNA/test.csv \
                       --localcores=8 \
                       --localmem=32