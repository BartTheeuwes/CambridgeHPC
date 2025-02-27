##!/bin/bash
eval "$(conda shell.bash hook)"
conda activate chip_peakcalling

bigWigToBedGraph /rds/project/rds-SDzz0CATGms/users/bt392/07_scRNA_CutTag/01_bulk/data/X204SC22051679-Z01-F001/raw_data/Cebpa/Cebpa.bw /rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood/ChIP/test.bedgraph