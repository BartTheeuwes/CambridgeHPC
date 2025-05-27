finemo extract-regions-bpnet-h5 -c EOd.profile_scores.h5 EOi.profile_scores.h5 \
    -o shap.npz
    
#finemo extract-regions-bpnet-h5 -c EOd.profile_scores.h5  -o shap.npz

finemo extract-regions-chrombpnet-h5 \
    -c ../PS_EOd/auxiliary/interpret_subsample/chrombpnet_nobias.profile_scores.h5 ../PS_EOi/auxiliary/interpret_subsample/chrombpnet_nobias.profile_scores.h5 \
    -o shap.npz
    
    
    
modisco motifs -s sequences.npz -a contributions.npz -n 20 -o modisco.h5 -w 500 -v

modisco report -i modisco.h5 -o modisco/ -m ~/miniconda3/envs/chrombpnet3/lib/python3.8/site-packages/chrombpnet/data/motifs.meme.txt

finemo call-hits -r shap.npz -m modisco.h5 -o finemo/