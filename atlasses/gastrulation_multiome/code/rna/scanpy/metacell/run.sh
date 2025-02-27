#!/nfs/research1/stegle/system/linuxbrew2/bin/zsh

################
## Define I/O ##
################

# if [ "$HOSTNAME" = ricard ]; then
# 	output_folder="/Users/ricard/data/gastrulation_multiome_10x/results/rna/metacells"
# else
# 	echo "Computer not recognised"; exit
# fi
output_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/results/rna/metacells/rna_atac"

samples=( "E7.5_rep1" "E7.5_rep2" "E8.0_rep1" "E8.0_rep2" "E8.5_rep1" "E8.5_rep2" )
# samples=( "E8.5_rep2" )

number_metacells=( 1000 2500 )
n_iter=50

for i in "${number_metacells[@]}"; do
	text_outfile="$output_folder/cell2metacell_${i}metacells.txt"
	anndata_outfile="$output_folder/anndata_${i}metacells.h5"

	cmd="python run_metacell.py --text_outfile $text_outfile --anndata_outfile $anndata_outfile --samples $samples --number_metacells $i --n_iter $n_iter"
	echo $cmd

	job 40 4 $cmd
	# eval $cmd

done

