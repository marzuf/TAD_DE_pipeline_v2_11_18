startTime <- Sys.time()
cat(paste0("> Rscript create_sameStrand_sortNoDup.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_sameStrand_sortNoDup.R TCGAcrc_msi_mss

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> add as.character() in apply !!! use "gene1" etc. instead of 1 index in apply 

caller <- "TopDom"

args <- commandArgs(trailingOnly = TRUE)
curr_dataset <- args[1]

stopifnot(length(args) == 1)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_SAME_STRAND_SORTNODUP", curr_dataset)
system(paste0("mkdir -p ", outFold))

dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

script0_name <- "0_prepGeneData"
pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
# 1/2-SBSRNA4         A2M       A2ML1       A2MP1      A4GALT       A4GNT 
# "100533182"         "2"    "144568"         "3"     "53947"     "51146" 

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
entrezDT$entrezID <- as.character(entrezDT$entrezID)

stopifnot(pipeline_geneList %in% entrezDT$entrezID)

geneDT <- entrezDT[entrezDT$entrezID %in% pipeline_geneList, c("entrezID", "strand")]
geneDT$entrezID <- as.character(geneDT$entrezID)
geneDT <- geneDT[order(geneDT$entrezID),]

all_gene_pairs <-  as.data.frame(t(combn(geneDT$entrezID, m=2)))
colnames(all_gene_pairs) <- c("gene1", "gene2")
all_gene_pairs$gene1 <- as.character(all_gene_pairs$gene1)
all_gene_pairs$gene2 <- as.character(all_gene_pairs$gene2)
stopifnot(all_gene_pairs$gene1 < all_gene_pairs$gene2)

all_strand_pairs <- merge(all_gene_pairs, geneDT, by.x="gene1", by.y="entrezID")
colnames(all_strand_pairs)[colnames(all_strand_pairs) == "strand"] <- "strand_gene1"

all_strand_pairs <- merge(all_strand_pairs, geneDT, by.x="gene2", by.y="entrezID")
colnames(all_strand_pairs)[colnames(all_strand_pairs) == "strand"] <- "strand_gene2"

all_strand_pairs$gene1 <- as.character(all_strand_pairs$gene1)
all_strand_pairs$gene2 <- as.character(all_strand_pairs$gene2)
stopifnot(all_strand_pairs$gene1 < all_strand_pairs$gene2)
head(all_strand_pairs)

all_strand_pairs$sameStrand <- as.numeric(all_strand_pairs$strand_gene1 == all_strand_pairs$strand_gene2)

all_strand_pairs$strand_gene1 <- all_strand_pairs$strand_gene2 <- NULL

head(all_strand_pairs)

outFile <- file.path(outFold, "all_strand_pairs.Rdata")
save(all_strand_pairs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))