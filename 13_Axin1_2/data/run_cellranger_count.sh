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
refdir=/rds/project/rds-SDzz0CATGms/references/10x/

# ====================================================================================





# =========================== EDIT THESE BELOW ===========================

# reference genome to use (select just one):
# ------------------------------------------
#ref=refdata-cellranger-hg19-1.2.0
#ref=refdata-cellranger-hg19_and_mm10-1.2.0  
ref=refdata-cellranger-mm10-1.2.0


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=/rds/project/rds-SDzz0CATGms/users/bt392/13_Axin1_2/data/
# libraries and respective expected sizes
# ---------------------------------------
library[LibA_DKOAxin1Axin2_E6_1_IGO_14590_1]=10000
#library[LibB_WT_E6_1_IGO_14590_2]=10000
# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/13_Axin1_2/data/cellranger_count.sh

# ========================================================================

echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
