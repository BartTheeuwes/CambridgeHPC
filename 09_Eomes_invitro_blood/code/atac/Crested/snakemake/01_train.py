# Train the crested model

import crested
import argparse
import os

# Check GPU availability
import tensorflow as tf
print("Num GPUs Available:", len(tf.config.list_physical_devices('GPU')))

# Define predefined fold configurations
FOLD_CONFIGS = {
    "fold0": {"val": 0.1, "test": 0.1, "seed": 42},
    "fold1": {"val": 0.12, "test": 0.09, "seed": 123},
    "fold2": {"val": 0.08, "test": 0.12, "seed": 2024},
    "fold3": {"val": 0.13, "test": 0.08, "seed": 99},
    "fold4": {"val": 0.09, "test": 0.13, "seed": 7},
}

parser = argparse.ArgumentParser(description="Import bigWigs data using crested.")
parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input folder containing bigwigs")
parser.add_argument("-p", "--peaks", type=str, required=True, help="Path to the peak file")
parser.add_argument("-o", "--output", type=str, required=True, help="Path to output folder")
parser.add_argument("-w", "--width", type=int, required=True, help="Target region width")
parser.add_argument("-f", "--fold", type=str, choices=FOLD_CONFIGS.keys(), required=True, help="Fold selection (fold0 to fold4)")

args = parser.parse_args()

# Retrieve corresponding val, test, and seed values
fold_config = FOLD_CONFIGS[args.fold]
val, test, seed = fold_config["val"], fold_config["test"], fold_config["seed"]

os.makedirs(os.path.join(args.output, args.fold), exist_ok=True)
os.chdir(os.path.join(args.output, args.fold))

# Call crested function with parsed arguments
adata = crested.import_bigwigs(
    bigwigs_folder=args.input,
    regions_file=args.peaks,
    target_region_width=args.width,
    target="count"
)

adata.obs.index = adata.obs.index.str.replace("_WT-TileSize-10-normMethod-ReadsInTSS-ArchR", "", regex=False)
adata = adata[[1,3,5,7,9,11,13,15,17,19,21],:]

print(adata)

crested.pp.train_val_test_split(
    adata,
    strategy="region",
    val_size=val,
    test_size=test,
    shuffle=True,
    random_state=seed,
)


crested.pp.change_regions_width(
    adata, 2114
)  # change the adata width of the regions to 2114bp


crested.pp.normalize_peaks(adata)

# Save the final preprocessing results
adata.write_h5ad(os.path.join(args.output, args.fold, 'preprocessed_data.h5ad'))

# Set the genome
genome = crested.Genome(
    "/home/bt392/rds/rds-bg200-hphi-gottgens/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa", 
    "/home/bt392/rds/rds-bg200-hphi-gottgens/users/bt392/10_Eomes_invitro_gut/results/atac/chrombpnet/mm10.chrom.sizes"
)
crested.register_genome(
    genome
)  # Register the genome so that it can be used by the package


datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,  # lower this if you encounter OOM errors
    max_stochastic_shift=3,  # optional data augmentation to slightly reduce overfitting
    always_reverse_complement=True,  # default True. Will double the effective size of the training dataset.
)


model_architecture = crested.tl.zoo.chrombpnet(seq_len=2114, num_classes=11)


# Load the default configuration for training a regression model
from crested.tl import default_configs, TaskConfig

config = default_configs("peak_regression")
print(config)


# setup the trainer
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=config,
    project_name = args.output,  # change to your liking
    run_name = args.fold,
    logger= None,  # or 'wandb', 'tensorboard'
)

# train the model
trainer.fit(epochs=25)

trainer.test()

# load the latest model
import keras
import glob
import re
# Get all .keras files in the specified directory
model_files = glob.glob(os.path.join(args.output, args.fold, "checkpoints/*.keras"))

# Extract numbers from filenames
def extract_number(filename):
    match = re.search(r"(\d+)\.keras", filename)
    return int(match.group(1)) if match else -1  # Return -1 if no number is found

# Find the latest model
latest_model_file = max(model_files, key=extract_number) if model_files else None

model = keras.models.load_model(
    os.path.join(args.output, args.fold, 'checkpoints', latest_model_file), compile=False
)

predictions = crested.tl.predict(adata, model)
adata.layers["predictions"] = predictions.T

import matplotlib.pyplot as plt


crested.pl.heatmap.correlations_self(
    adata, title="Self Correlation Heatmap", x_label_rotation=90, width=6, height=6
)
plt.savefig(args.output + '/' + args.fold + '/self_correlations_plot.pdf', format='pdf')
plt.close()

crested.pl.heatmap.correlations_predictions(
    adata,
    split="train",
    title="Correlations between Groundtruths and Train Predictions",
    x_label_rotation=90,
    width=5,
    height=5,
    log_transform=True,
    vmax=1,
    vmin=-0.15,
)
plt.savefig(args.output + '/' + args.fold + '/train_correlations_predictions_plot.pdf', format='pdf')
plt.close()

crested.pl.heatmap.correlations_predictions(
    adata,
    split="val",
    title="Correlations between Groundtruths and Val Predictions",
    x_label_rotation=90,
    width=5,
    height=5,
    log_transform=True,
    vmax=1,
    vmin=-0.15,
)
plt.savefig(args.output + '/' + args.fold + '/val_correlations_predictions_plot.pdf', format='pdf')
plt.close()

crested.pl.heatmap.correlations_predictions(
    adata,
    split="test",
    title="Correlations between Groundtruths and Test Predictions",
    x_label_rotation=90,
    width=5,
    height=5,
    log_transform=True,
    vmax=1,
    vmin=-0.15,
)
plt.savefig(args.output + '/' + args.fold + '/test_correlations_predictions_plot.pdf', format='pdf')
plt.close()

# script completion file
with open(os.path.join(args.output, args.fold, "script_completed.txt"), "w") as f:
    f.write("Script finished successfully.\n")







