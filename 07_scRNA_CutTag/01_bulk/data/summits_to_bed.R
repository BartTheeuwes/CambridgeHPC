library(tidyr)
library(data.table)

TFs = c('Cebpa', 'Lyl1', 'H3AcK27', 'Pu1')
for(i in TFs){

	TF_i = fread(sprintf('/rds/project/rds-SDzz0CATGms/users/bt392/07_scRNA_CutTag/01_bulk/data/X204SC22051679-Z01-F001/raw_data/%s/macs2/%s_summits.bed', i, i))[,c(1,2,3)] %>%
		setnames(c('chr', 'start', 'end'))
	TF_bed = TF_i %>% .[,start:=start-250] %>% .[,end:=end+250]
	fwrite(TF_bed, sprintf('/rds/project/rds-SDzz0CATGms/users/bt392/07_scRNA_CutTag/01_bulk/data/X204SC22051679-Z01-F001/raw_data/%s/macs2/%s_500bp.bed', i, i), sep = "\t", quote = FALSE, col.names = FALSE)
	
}