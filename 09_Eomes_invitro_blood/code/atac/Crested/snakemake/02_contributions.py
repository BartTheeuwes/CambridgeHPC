# load the latest model
import keras
import glob
import re

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
parser.add_argument("-ct", "--celltype", type=str, required=True, help="cell type to profile")


args = parser.parse_args()

os.makedirs(os.path.join(args.output, args.runname, 'contributions'), exist_ok=True)

# Get all .keras files in the specified directory
model_files = glob.glob(os.path.join(args.output, args.fold, "checkpoints/*.keras"))

# Extract numbers from filenames
def extract_number(filename):
    match = re.search(r"(\d+)\.keras", filename)
    return int(match.group(1)) if match else -1  # Return -1 if no number is found

# Find the latest model
latest_model_file = max(model_files, key=extract_number) if model_files else None

model = keras.models.load_model(
    latest_model_file, compile=False
)

# Call crested function with parsed arguments
adata = crested.import_bigwigs(
    bigwigs_folder=args.input,
    regions_file=args.peaks,
    target_region_width=args.width,
    target="count"
)

adata.obs.index = adata.obs.index.str.replace("_WT-TileSize-10-normMethod-ReadsInTSS-ArchR", "", regex=False)
adata = adata[[1,3,5,7,9,11,13,15,17,19,21],:]

# Set the genome
genome = crested.Genome(
    "/home/bt392/rds/rds-bg200-hphi-gottgens/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa", 
    "/home/bt392/rds/rds-bg200-hphi-gottgens/users/bt392/10_Eomes_invitro_gut/results/atac/chrombpnet/mm10.chrom.sizes"
)
crested.register_genome(
    genome
)  # Register the genome so that it can be used by the package


# extract index number for cell type
index_number = adata.obs.index.get_loc(args.celltype)

crested.tl.contribution_scores(   
    input=adata,
    target_idx=index_number,  # We calculate for all classes
    model=model,
    output_dir=os.path.join(args.output, args.fold, 'contributions'),
    method="integrated_grad"
)
