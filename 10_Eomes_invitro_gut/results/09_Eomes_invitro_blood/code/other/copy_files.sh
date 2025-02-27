
#################
## First batch ##
#################

# eomes_2_L001: rv_eo_deg_day3_5_dtag
# eomes_3_L001: rv_eo_deg_day3_5_control
# eomes_4_L001: rv_eo_deg_day4_dtag
# eomes_5_L001: rv_eo_deg_day4_control
# eomes_6_L001: rv_eo_deg_day4_5_dtag
# eomes_7_L001: rv_eo_deg_day4_5_control
# eomes_8_L001: eo_ko_d4_5

basedir="/bi/sequencing/Sample_5568_Multiome_eomes/Lane_7940_Multiome_eomes/"
samples=( "eomes_2_L001" "eomes_3_L001" "eomes_4_L001" "eomes_5_L001" "eomes_6_L001" "eomes_7_L001" "eomes_8_L001")
samples_alias=( "rv_eo_deg_day3_5_dtag" "rv_eo_deg_day3_5_control" "rv_eo_deg_day4_dtag" "rv_eo_deg_day4_control" "rv_eo_deg_day4_5_dtag" "rv_eo_deg_day4_5_control" "eo_ko_d4_5")

ln -s ${basedir}/Aligned/*gex_possorted_bam.bam .
ln -s ${basedir}/Unaligned/*gex_possorted_bam.bam.bai .
ln -s ${basedir}/Unaligned/*atac_fragments.tsv.gz .
ln -s ${basedir}/Aligned/*barcodes.tsv.gz .
ln -s ${basedir}/Unaligned/*features.tsv.gz .
ln -s ${basedir}/Aligned/*matrix.mtx.gz .
ln -s ${basedir}/Aligned/*web_summary.html .

rename "eomes_2_L001" "rv_eo_deg_day3_5_dtag" eomes_2_L001*
rename "eomes_3_L001" "rv_eo_deg_day3_5_control" eomes_3_L001*
rename "eomes_4_L001" "rv_eo_deg_day4_dtag" eomes_4_L001*
rename "eomes_5_L001" "rv_eo_deg_day4_control" eomes_5_L001*
rename "eomes_6_L001" "rv_eo_deg_day4_5_dtag" eomes_6_L001*
rename "eomes_7_L001" "rv_eo_deg_day4_5_control" eomes_7_L001*
rename "eomes_8_L001" "eo_ko_d4_5" eomes_8_L001*

for i in ${samples_alias[@]}; do
	mkdir -p $i/outs
	mv "${i}_gex_possorted_bam.bam" $i/outs
	mv "${i}_gex_possorted_bam.bam.bai" $i/outs
	mv "${i}_web_summary.html" $i/outs
	mv "${i}_atac_fragments.tsv.gz" $i/outs

	rename "${i}_" "" $i/outs/*

	mkdir $i/outs/unfiltered_feature_bc_matrix
	mv "${i}_unfiltered_barcodes.tsv.gz" $i/outs/unfiltered_feature_bc_matrix/
	mv "${i}_unfiltered_features.tsv.gz" $i/outs/unfiltered_feature_bc_matrix/
	mv "${i}_unfiltered_matrix.mtx.gz" $i/outs/unfiltered_feature_bc_matrix/
	rename "${i}_unfiltered_" "" $i/outs/unfiltered_feature_bc_matrix/*


	mkdir $i/outs/filtered_feature_bc_matrix
	mv "${i}_barcodes.tsv.gz" $i/outs/filtered_feature_bc_matrix/
	mv "${i}_features.tsv.gz" $i/outs/filtered_feature_bc_matrix/
	mv "${i}_matrix.mtx.gz" $i/outs/filtered_feature_bc_matrix/
	rename "${i}_" "" $i/outs/filtered_feature_bc_matrix/*
	
done


##################
## Second batch ##
##################

# (everything else) /bi/sequencing/Sample_5688/Lane_8167_SLX-20897_Eomes_multiome_batch2_GEX
# (fastqc atac) /bi/sequencing/Sample_5687_SLX-20898_Eomes_multiome_batch2_ATAC/Lane_8168_SLX-20898_Eomes_multiome_batch2_ATAC

basedir="/bi/sequencing/Sample_5688_SLX-20897_Eomes_multiome_batch2_GEX/Lane_8167_SLX-20897_Eomes_multiome_batch2_GEX"
samples=( "1A_Eo_DEG_G9_day3" "1B_Eo_DEG_G9_day3" "2_Eo_DEG_G9_day3_5_dTAG" "2_Eo_DEG_G9_day3_5_VC" "2_Eo_DEG_G9_day4_dTAG" "2_Eo_DEG_G9_day4_VC" "2_Eo_DEG_G9_day5_dTAG" "2_Eo_DEG_G9_day5_VC" )


# i="1A_Eo_DEG_G9_day3"
for i in ${samples[@]}; do
	echo "$i"
	mkdir -p $i/outs
	cd $i/outs

	ln -s "${basedir}/Aligned/${i}_L001_gex_possorted_bam.bam" "gex_possorted_bam.bam"
	ln -s "${basedir}/Unaligned/${i}_L001_gex_possorted_bam.bam.bai" "gex_possorted_bam.bam.bai"
	ln -s "${basedir}/Aligned/${i}_L001_web_summary.html" "web_summary.html"
	ln -s "${basedir}/Unaligned/${i}_L001_atac_fragments.tsv.gz" "atac_fragments.tsv.gz"
	# rename "${i}_L001_" "" *

	mkdir filtered_feature_bc_matrix
	cd filtered_feature_bc_matrix

	ln -s "${basedir}/Aligned/${i}_L001_barcodes.tsv.gz"  "barcodes.tsv.gz"
	ln -s "${basedir}/Unaligned/${i}_L001_features.tsv.gz"  "features.tsv.gz"
	ln -s "${basedir}/Aligned/${i}_L001_matrix.mtx.gz"  "matrix.mtx.gz"
	# rename "${i}_L001_" "" *

	cd ..
	mkdir unfiltered_feature_bc_matrix
	cd unfiltered_feature_bc_matrix

	ln -s "${basedir}/Aligned/${i}_L001_unfiltered_barcodes.tsv.gz" "barcodes.tsv.gz"
	ln -s "${basedir}/Unaligned/${i}_L001_unfiltered_features.tsv.gz" "features.tsv.gz"
	ln -s "${basedir}/Aligned/${i}_L001_unfiltered_matrix.mtx.gz" "matrix.mtx.gz"
	# rename "${i}_L001_" "" *

	cd ../../..
done

##########
## TEST ##
##########

for i in ${samples[@]}; do
	echo "$i"
	# mv $i/outs/unfiltered_feature_bc_matrix/matrix.tsv.gz $i/outs/unfiltered_feature_bc_matrix/matrix.mtx.gz
	ln -fs ${basedir}/Aligned/${i}_L001_unfiltered_matrix.mtx.gz $i/outs/unfiltered_feature_bc_matrix/matrix.mtx.gz
done


# Files needed:
# - atac_fragments.tsv.gz
# - atac_fragments.tsv.gz.tbi
# - barcodes.tsv.gz
# - features.tsv.gz
# - matrix.mtx.gz
# - web_summary.html