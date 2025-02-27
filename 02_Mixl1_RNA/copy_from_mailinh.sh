#!/bin/bash
#SBATCH -p skylake #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl3-cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time 12:00:00
#SBATCH --job-name copy

cp -r /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/PhD_MT07_Mixl1/SLX19260_Mixl1_SIGAA5 data/
cp -r /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/PhD_MT07_Mixl1/SLX19260_Mixl1_SIGAB5 data/
cp -r /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/PhD_MT07_Mixl1/SLX19260_Mixl1_SIGAC5 data/
cp -r /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/PhD_MT07_Mixl1/SLX19260_Mixl1_SIGAD5 data/
cp -r /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/PhD_MT07_Mixl1/SLX19260_Mixl1_SIGAE5 data/
cp -r /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/PhD_MT07_Mixl1/SLX19260_Mixl1_SIGAF5 data/