# TO DO:
# Script to mean contributions of 5 folds:
# https://github.com/kundajelab/neuro-variants/blob/main/preprocess/4_shap_peaks/trevino_2021/1_get_mean_shap_scores.py

Or could use finemo:???

https://github.com/austintwang/finemo_gpu
finemo extract-regions-chrombpnet-h5
Extract sequences and contributions from ChromBPNet H5 files.

Usage: finemo extract-regions-chrombpnet-h5 -c <h5s> -o <out_path> [-w <region_width>]

-c/--h5s: One or more H5 files of contribution scores, with paths delimited by whitespace. Scores are averaged across files.

# Could this output be used as an input for 08_TFmodisco.sh??
# output = average of all contribution bigwigs, and saved as compressed .npz
# can the npz be used as direct input, or just need to decompress back to h5 form TFmodisco?
# npz can be used but might need to remove the one-hot first?