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

indir=/rds/project/rds-SDzz0CATGms/users/bt392/06_Runx1_RNA/data/SLX21184/


# ================== PROBABLY THERE WILL BE NO NEED TO CHANGE THESE ==================

declare -A library
refdir=/rds/project/rds-SDzz0CATGms/references/10x/

# ====================================================================================





# =========================== EDIT THESE BELOW ===========================

# reference genome to use (select just one):
# ------------------------------------------
#ref=refdata-cellranger-hg19-1.2.0
#ref=refdata-cellranger-hg19_and_mm10-1.2.0  
ref=refdata-gex-mm10-2020-A


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=$indir'genotyping/'
# libraries and respective expected sizes
# ---------------------------------------
library[SLX-21184_SITTA7_H5NKNDMXY]=13784
library[SLX-21184_SITTB7_H5NKNDMXY]=946
library[SLX-21184_SITTC7_H5NKNDMXY]=13784
library[SLX-21184_SITTD7_H5NKNDMXY]=946
library[SLX-21184_SITTE12_H5NKNDMXY]=13784
library[SLX-21184_SITTH10_H5NKNDMXY]=946


# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=$fastqdir'cellranger_count.sh'
# ========================================================================




echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
