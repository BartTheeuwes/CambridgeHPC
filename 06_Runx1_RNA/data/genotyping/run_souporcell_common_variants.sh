#!/bin/bash
souporcell="/rds/project/rds-SDzz0CATGms/users/bt392/06_Runx1_RNA/data/genotyping/souporcell_common_variants.sh"

declare -A library
library[SITTA7]=
library[SITTB7]=
library[SITTC7]=
library[SITTD7]=
library[SITTE12]=
library[SITTH10]=

for i in "${!library[@]}"
do
  echo ${i}
  sbatch ${souporcell} ${i}
done