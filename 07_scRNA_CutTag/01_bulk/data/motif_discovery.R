suppressPackageStartupMessages({
    library(rGADEM)
    library(data.table)
    library(tidyr)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

TF = 'Cebpa'

TF.gr = fread("/rds/project/rds-SDzz0CATGms/users/bt392/07_scRNA_CutTag/01_bulk/data/X204SC22051679-Z01-F001/raw_data/Cebpa/macs2/NA_summits.bed") %>%
    .[,c('V1', 'V2', 'V3')] %>% 
    setnames(c('chr', 'start', 'end')) %>%
    makeGRangesFromDataFrame %>% resize(., width = 100, fix='center')

TF.seq = getSeq(Mmusculus, TF.gr)

GADEM = GADEM(TF.seq, 
              seed=1234,
              genome=Mmusculus,
              nmotifs=10)

saveRDS(GADEM, "/rds/project/rds-SDzz0CATGms/users/bt392/07_scRNA_CutTag/01_bulk/data/X204SC22051679-Z01-F001/raw_data/Cebpa/macs2/GADEM.rds"))