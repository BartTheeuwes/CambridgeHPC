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




# ================== PROBABLY THERE WILL BE NO NEED TO CHANGE THESE ==================

declare -A library
refdir=/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/

# ====================================================================================





# =========================== EDIT THESE BELOW ===========================

# reference genome to use (select just one):
# ------------------------------------------
#ref=refdata-cellranger-hg19-1.2.0
#ref=refdata-cellranger-hg19_and_mm10-1.2.0  
ref=mm10_tomato


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/raw_data/
# libraries and respective expected sizes
# ---------------------------------------
library[SIGAA5]=4319
library[SIGAB5]=4319
library[SIGAC5]=4319
library[SIGAD5]=4139
library[SIGAE5]=4319
library[SIGAF5]=4319

# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/raw_data/cellranger_count.sh
# ========================================================================




echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
