#!/bin/bash

conda_mv="/rds/project/rds-SDzz0CATGms/users/bt392/Conda_mv.sh"
conda_mv_file="/rds/project/rds-SDzz0CATGms/users/bt392/Conda_mv_file.sh"

mv /home/bt392/miniconda3/*.txt /rds/project/rds-SDzz0CATGms/users/bt392/miniconda3/

cd /home/bt392/miniconda3/
mainfolder='sbin' #'ls -d */'
for folder in $mainfolder
do
  mkdir /rds/project/rds-SDzz0CATGms/users/bt392/miniconda3/$folder
  # move files
  cd /home/bt392/miniconda3/$folder
  sbatch $conda_mv_file $mainfolder
  # move subfolders 
  subfolders='ls -d */'
  for subfolder in $subfolders
  do
    echo $mainfolder $subfolder
  done
done