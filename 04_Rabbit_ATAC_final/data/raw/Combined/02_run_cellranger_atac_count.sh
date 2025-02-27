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

# RUN THIS MANUALLY IN TERMINAL!
export PATH=/rds/project/rds-SDzz0CATGms/users/bt392/software/cellranger-atac-2.0.0/cellranger-atac-2.0.0:$PATH

# ================== PROBABLY THERE WILL BE NO NEED TO CHANGE THESE ==================

declare -A library
refdir=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/genome/10x_rabbit_ref/

# ====================================================================================





# =========================== EDIT THESE BELOW ===========================

# reference genome to use (select just one):
# ------------------------------------------
#ref=refdata-cellranger-hg19-1.2.0
#ref=refdata-cellranger-hg19_and_mm10-1.2.0  
ref=OryCun2


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/data/raw/SLX18911/
# libraries and respective expected sizes
# ---------------------------------------
#library[SLX-20795_SITTG10_HKTG2DRXY]=4465
library[SLX18911_SINAA3,SLX19158_]=13784
library[SLX18911_SINAB3,SLX19158_]=946
library[SLX18911_SINAC3,SLX19158_]=13784
library[SLX18911_SINAD3,SLX19158_]=946
library[SLX18911_SINAE3,SLX19158_]=13784
library[SLX18911_SINAF3,SLX19158_]=946
library[SLX18911_SINAG3,SLX19158_]=13784
library[SLX18911_SINAH3,SLX19158_]=946

# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/data/raw/SLX18911/cellranger_atac_count.sh
# ========================================================================




echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
