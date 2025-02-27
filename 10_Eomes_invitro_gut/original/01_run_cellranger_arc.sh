#!/bin/bash
# -------------------------------------------------------
#title           :run_cellranger_count.sh
#author          :hpb29
#date            :20180404
#version         :0.5    
#description     :Controls the submition of 10x libraries
#                 via SLURM to 'cellranger count'
#usage           :Edit variables below accordingly and run
# --------------------------------------------------------

# libraries and respective expected sizes
# ---------------------------------------
declare -A library
#library[RBG43471]=13784
library[RBG43472]=13784
library[RBG43473]=13784
library[RBG43474]=13784
library[RBG43475]=13784
library[RBG43476]=13784
library[RBG43477]=13784
library[RBG43478]=13784

CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/10_Eomes_invitro_gut/original/01_cellranger_arc.sh

echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i}
done
