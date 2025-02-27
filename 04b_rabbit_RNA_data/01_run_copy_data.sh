#!/bin/bash
mkdir B1
mkdir B2
mkdir B3
mkdir B4

copy_data=/rds/project/rds-SDzz0CATGms/users/bt392/04b_rabbit_RNA_data/01_copy_data.sh
search_dir=/rfs/project/rfs-jlyrZBLdWBU/Sequencing_raw_data/Gottgens/SLX18995/
for entry in "$search_dir"/*.tar
do
  sbatch ${copy_data} $entry B1
done

search_dir=/rfs/project/rfs-jlyrZBLdWBU/Sequencing_raw_data/Gottgens/SLX18995_B2/
for entry in "$search_dir"/*.tar
do
  sbatch ${copy_data} $entry B2
done

search_dir=/rfs/project/rfs-jlyrZBLdWBU/Sequencing_raw_data/Gottgens/SLX18995_B3/
for entry in "$search_dir"/*.tar
do
  sbatch ${copy_data} $entry B3
done

search_dir=/rfs/project/rfs-jlyrZBLdWBU/Sequencing_raw_data/Gottgens/SLX18995_B4/
for entry in "$search_dir"/*.tar
do
  sbatch ${copy_data} $entry B4
done


