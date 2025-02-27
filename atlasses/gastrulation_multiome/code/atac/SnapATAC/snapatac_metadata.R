library(SnapATAC)
library(purrr)
library(data.table)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))



snap <- readRDS(snapio$rds_file)
snap

#sample_metadata <- fread(io$metadata)

snapmeta <- setDT(snap@metaData)

# add umap coordinates
snapmeta <- cbind(snapmeta, as.data.table(snap@umap)) 
  

# add chromVar data
snapmeta <- cbind(snapmeta, as.data.table(snap@mmat))

# re-format names
make_names <- function(x){
  make.unique(x, sep = "_") %>%
    gsub("-|\\.|\\::|\\(|\\)", "_", .) %>% 
    gsub("var_2_", "var2", .)
}

newnames <- make_names(colnames(snapmeta))
dt=data.table(colnames(snapmeta), newnames)
dt[V1 %like% "\\::" | V1 %like% "\\("]

colnames(snapmeta) <- newnames
snapmeta

fwrite(snapmeta, snapio$metadata, sep = "\t", na = "NA", quote = FALSE)
