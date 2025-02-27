library(rhdf5)
ArrowFile <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/processed/atac/archR_integrated/ArrowFiles/E8.25_PijuanSala.arrow"

fid <- H5Fopen(ArrowFile)
h5ls(ArrowFile)
h5read(fid, "Metadata")
h5read(fid, "Metadata/CellNames")
h5read(fid, "Metadata/Sample")
h5write("E8.25_PijuanSala", file=fid, name="Metadata/Sample")

h5closeAll()
