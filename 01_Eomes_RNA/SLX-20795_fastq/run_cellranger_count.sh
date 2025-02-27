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
ref=mm10_tomato


# path to dir where 10x fastq files are
# --------------------------------------------------------------
fastqdir=/rds/project/rds-SDzz0CATGms/users/bt392/phd_01_Eomes_RNA/SLX-20795_fastq/
# libraries and respective expected sizes
# ---------------------------------------
#library[SLX-20795_SITTG10_HKTG2DRXY]=4465
library[SLX-20795_SITTH10_HKTG2DRXY]=13784
library[SLX-20795_SITTA11_HKTG2DRXY]=946
library[SLX-20795_SITTB11_HKTG2DRXY]=8573
library[SLX-20795_SITTF11_HKTG2DRXY]=4561
library[SLX-20795_SITTG11_HKTG2DRXY]=4745
library[SLX-20795_SITTH11_HKTG2DRXY]=1513
library[SLX-20795_SITTA12_HKTG2DRXY]=3463

# (...)

# path to cellranger_count.sh file
CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/phd_01_Eomes_RNA/SLX-20795_fastq/cellranger_count.sh
# ========================================================================




echo "Resulting analysis pipestance folders will be output at:"
echo " `pwd`"

for i in "${!library[@]}"
do
  sbatch ${CELLRANGERCTRL} ${i} ${refdir}${ref} ${fastqdir} ${library["$i"]}  
done
