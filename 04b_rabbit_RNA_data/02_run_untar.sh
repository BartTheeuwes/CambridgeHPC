#!/bin/bash
untar=/rds/project/rds-SDzz0CATGms/users/bt392/04b_rabbit_RNA_data/02_untar.sh
search_dir=/rds/project/rds-SDzz0CATGms/users/bt392/04b_rabbit_RNA_data/B1/
for entry in "$search_dir"/*.tar
do
 sbatch ${untar} $entry $search_dir
done

search_dir=/rds/project/rds-SDzz0CATGms/users/bt392/04b_rabbit_RNA_data/B2/
for entry in "$search_dir"/*.tar
do
 sbatch ${untar} $entry $search_dir
done

search_dir=/rds/project/rds-SDzz0CATGms/users/bt392/04b_rabbit_RNA_data/B3/
for entry in "$search_dir"/*.tar
do
 sbatch ${untar} $entry $search_dir
done

search_dir=/rds/project/rds-SDzz0CATGms/users/bt392/04b_rabbit_RNA_data/B4/
for entry in "$search_dir"/*.tar
do
 sbatch ${untar} $entry $search_dir
done