# RUN THIS MANUALLY IN TERMINAL!
export PATH=/rds/project/rds-SDzz0CATGms/users/bt392/software/cellranger-atac-2.0.0/cellranger-atac-2.0.0:$PATH

# ================== PROBABLY THERE WILL BE NO NEED TO CHANGE THESE ==================

declare -A library
refdir=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/genome/10x_rabbit_ref/
ref=OryCun2
fastqdir=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/data/raw/Combined/
CELLRANGERCTRL=/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/data/raw/Combined/cellranger_atac_count.sh

sbatch ${CELLRANGERCTRL} BGRGP1 ${refdir}${ref} ${fastqdir} SLX18911_SINAA3,SLX19158_SINAA3 
sbatch ${CELLRANGERCTRL} BGRGP2 ${refdir}${ref} ${fastqdir} SLX18911_SINAB3,SLX19158_SINAB3
sbatch ${CELLRANGERCTRL} BGRGP3 ${refdir}${ref} ${fastqdir} SLX18911_SINAC3,SLX19158_SINAC3  
sbatch ${CELLRANGERCTRL} BGRGP4 ${refdir}${ref} ${fastqdir} SLX18911_SINAD3,SLX19158_SINAD3 
sbatch ${CELLRANGERCTRL} BGRGP5 ${refdir}${ref} ${fastqdir} SLX18911_SINAE3,SLX19158_SINAE3  
sbatch ${CELLRANGERCTRL} BGRGP6 ${refdir}${ref} ${fastqdir} SLX18911_SINAF3,SLX19158_SINAF3
sbatch ${CELLRANGERCTRL} BGRGP7 ${refdir}${ref} ${fastqdir} SLX18911_SINAG3,SLX19158_SINAG3
sbatch ${CELLRANGERCTRL} BGRGP8 ${refdir}${ref} ${fastqdir} SLX18911_SINAH3,SLX19158_SINAH3