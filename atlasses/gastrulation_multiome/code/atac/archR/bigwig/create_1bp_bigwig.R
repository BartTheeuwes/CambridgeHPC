suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(parallel))

here::i_am("atac/archR/bigwig/create_1bp_bigwigs.ipynb")

source(here::here("settings.R"))
source(here::here("utils.R"))

ArchRProject = loadArchRProject(io$archR.directory)

meta = as.data.table(ArchRProject@cellColData, keep.rownames = T)

celltypes = unique(ArchRProject@cellColData$celltype.mapped)
celltypes = c('Endothelium', 'Primitive_Streak', 'Erythroid_1', 'Erythroid_3', 'Allantois', 'Nascent_Mesoderm')

samples = list.files('/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/data/original/')
samples = samples[!grepl('CRISPR', samples)]

chrs = c('chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX','chrY')

lapply(celltypes, function(x){
    message(paste0('Celltype: ', x))

    full = data.table(chr = NA, position = NA, reads = NA)

    for(i in samples){
        tmp_meta = meta[celltype.mapped == x & sample == i] %>% 
            .[,cell := str_split(rn, '#') %>% map_chr(2)]
        
        if(nrow(tmp_meta) >= 100){
            message(paste0('reading sample: ', i))
            a = Sys.time()
            file = sprintf('/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/data/original/%s/atac_fragments.tsv.gz', i)
            fragment = suppressWarnings(fread(file,
                            tmpdir = '/rds/project/rds-SDzz0CATGms/users/bt392/software/tmp'))
            b = Sys.time()
            message(paste0('File read time: ', b-a))

            fragment = fragment %>% 
                setnames(c('chr', 'start', 'end', 'cell', 'reads')) %>%
                # .[chr == 'chr1'] %>%
                .[cell %in% tmp_meta[,cell]] %>% 
                .[chr %in% chrs]

            tmp = mclapply(1:nrow(fragment), function(x){
                read = fragment[x,]
                tmp = data.table(chr = read$chr,
                                 position = seq(read$start, read$end, 1),
                                 reads = 1)
                return(tmp)
            }, mc.cores = 36) %>% rbindlist()

            full = rbind(full, tmp) %>% 
                .[, reads := sum(reads), by = c('chr', 'position')] %>% unique(by = c('chr', 'position'))
        }
    }
    
    full = full[-1] # Remove NA value at the very beginning
    
    message('Summing over ranges')
    full_summed = mclapply(chrs, function(chromosome){
        tmp = full[chr == chromosome][order(position)]

        # Collapse consecutive rows where the 'group' column is the same
        collapsed_dt <- tmp[, .(
          start = min(position),       # Starting id of the collapsed group
          end = max(position)          # Ending id of the collapsed group
        ), by = .(reads, rleid(reads))] %>% 
            .[,chr := chromosome] %>%
            .[,.(chr, start, end, reads)] %>% 
            setnames(c('chr', 'start', 'end', 'score'))

        return(collapsed_dt)
    }, mc.cores = length(chrs)) %>% rbindlist() 
              
    message('Making bigwigs')
    tmp_bw = makeGRangesFromDataFrame(full_summed, keep.extra.columns = T)
    tmp_bw@seqinfo = ArchRProject@geneAnnotation$genes@seqinfo

    bw_zeros = GenomicRanges::gaps(tmp_bw)
    bw_zeros = bw_zeros[bw_zeros@strand == '*']
    bw_zeros$score = 0

    tmp = rbind(as.data.table(tmp_bw), as.data.table(bw_zeros)) %>% 
        .[order(seqnames, start)]
    
    bw = makeGRangesFromDataFrame(tmp, keep.extra.columns = T)
    bw$score = bw$score / sum(meta[celltype.mapped == x, nFrags])
    bw@seqinfo = ArchRProject@geneAnnotation$genes@seqinfo
                         
    output_bw_file = sprintf('/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/data/processed/atac/archR/GroupBigWigs/test/%s.bw', x)
    rtracklayer::export(bw, output_bw_file, format = "BigWig")
})
