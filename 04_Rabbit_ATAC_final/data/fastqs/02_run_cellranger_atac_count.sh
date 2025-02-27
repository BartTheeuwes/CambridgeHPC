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
refdir=/rds/project/rds-SDzz0CATGms/users/mlnt2/rabbit_reference/

# ====================================================================================





# =========================== EDIT THESE BELOW ===========================

# reference genome to use (select just one):
# ------------------------------------------
#ref=refdata-cellranger-hg19-1.2.0
#ref=refdata-cellranger-hg19_and_mm10-1.2.0  
ref=ocun_atac


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/data/raw/SLX18911/
# libraries and respective expected sizes
# ---------------------------------------
#library[SLX-20795_SITTG10_HKTG2DRXY]=4465
library[SINAA3]=13784
library[SINAB3]=946
library[SINAC3]=13784
library[SINAD3]=946
library[SINAE3]=13784
library[SINAF3]=946
library[SINAG3]=13784
library[SINAH3]=946

# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/SLX-21143_fastq/cellranger_atac_count.sh
# ========================================================================

export PATH=/opt/cellranger-atac-2.0.0:$PATH


echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
