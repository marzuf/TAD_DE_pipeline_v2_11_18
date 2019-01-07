# Rscript get_tad_symbols.R chr19_TAD28

setDir=""

gene2tad_dt <- read.delim(file.path(setDir, 
                                    "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt"),
                          col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE, header=F)

gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)


entrez2symb_dt <- read.delim(file.path(setDir,
                                      "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)

library(dplyr)

mergeDT <- left_join(gene2tad_dt, entrez2symb_dt, by ="entrezID")


args <- commandArgs(trailingOnly = T)
curr_tad = args[1]


tad_dt <- mergeDT[mergeDT$region == curr_tad,c("region", "entrezID", "symbol")]


write.table(tad_dt, col.names=T, row.names=F, quote=F, sep="\t", file="")
