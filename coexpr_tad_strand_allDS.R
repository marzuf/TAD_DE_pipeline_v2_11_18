# if promoter sharing as a source of coexpression,
# divergently transcriped genes should have higher coexpr.

library(foreach)
library(doMC)

setDir <- "/media/electron"
registerDoMC(2)
setDir <- ""
registerDoMC(40)

source("analysis_utils.R")

plotType <- "png"
myHeight <- 400
myWidth <- 600

outFold <- "COEXPR_TAD_STRAND_allDS"
system(paste0("mkdir -p ", outFold))

curr_dataset <- "TCGAcrc_msi_mss"
myTit <- paste0("All datasets", " - Coexpression: TAD and strand ")

all_datasets <- list.files(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"))
nDS <- length(all_datasets)

all_ds_values <- foreach(curr_dataset = all_datasets, .combine='rbind') %dopar% {
  
  cat("... start: ", curr_dataset, "\n")
  
  coexprFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUC_COEXPRDIST_SORTNODUP", curr_dataset, "allData_dt.Rdata")
  
  coexprDT <- eval(parse(text=load(coexprFile)))
  coexprDT$gene1 <- as.character(coexprDT$gene1)
  coexprDT$gene2 <- as.character(coexprDT$gene2)
  
  entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
  entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
  entrezDT$entrezID <- as.character(entrezDT$entrezID)
  
  stopifnot(coexprDT$gene1 %in% entrezDT$entrezID)
  stopifnot(coexprDT$gene2 %in% entrezDT$entrezID)
  
  
  coexprDT <- merge(coexprDT, entrezDT[,c("entrezID", "strand")], by.x = "gene1", by.y="entrezID")
  colnames(coexprDT)[colnames(coexprDT) == "strand"] <- "strand_gene1"
  
  coexprDT <- merge(coexprDT, entrezDT[,c("entrezID", "strand")], by.x = "gene2", by.y="entrezID")
  colnames(coexprDT)[colnames(coexprDT) == "strand"] <- "strand_gene2"
  
  coexprDT$strands <- ifelse(coexprDT$strand_gene1 == coexprDT$strand_gene2, "sameStrand", "diffStrand")
  
  coexprDT$sameTAD <- as.character(coexprDT$sameTAD)
  coexprDT$sameTAD[coexprDT$sameTAD == 1] <- "sameTAD"
  coexprDT$sameTAD[coexprDT$sameTAD == 0] <- "diffTAD"
  coexprDT$sameTAD <- factor(coexprDT$sameTAD, levels = c("sameTAD", "diffTAD"))
  
  coexprDT$strands <- factor(coexprDT$strands, levels = c("diffStrand", "sameStrand"))
  
  # outFile <- file.path(outFold, paste0(curr_dataset, "_boxplot_coexpr_TAD_strand.svg"))
  # do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  # boxplot(coexpr ~ strands + sameTAD, data = coexprDT, main = myTit,
  #         ylab = "Gene coexpression (PCC)")
  # 
  # foo <- dev.off()
  # cat(paste0("... written: ", outFile, "\n"))
  
  # table(coexprDT$strands, coexprDT$sameTAD)
  
  
  # outFile <- file.path(outFold, paste0(curr_dataset, "_multiline_density_coexpr_TAD_strand.svg"))
  # do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  # plot_multiDens(list(
  #   divergent_sameTAD = coexprDT$coexpr[coexprDT$sameTAD == "sameTAD" & coexprDT$strands == "diffStrand"],
  #   same_sameTAD = coexprDT$coexpr[coexprDT$sameTAD == "sameTAD" & coexprDT$strands == "sameStrand"],
  #   divergent_diffTAD = coexprDT$coexpr[coexprDT$sameTAD == "diffTAD" & coexprDT$strands == "diffStrand"],
  #   same_diffTAD = coexprDT$coexpr[coexprDT$sameTAD == "diffTAD" & coexprDT$strands == "sameStrand"]
  # ), my_xlab =  "Gene coexpression (PCC)", plotTit  = myTit
  # )
  # foo <- dev.off()
  # cat(paste0("... written: ", outFile, "\n"))
  
  coexprDT[, c("coexpr", "sameTAD", "strands")]
  
}

outFile <- file.path(outFold, "all_ds_values.Rdata")
save(all_ds_values, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# load("COEXPR_TAD_STRAND_allDS/all_ds_values.Rdata")

outFile <- file.path(outFold, paste0("all_datasets", "_boxplot_coexpr_TAD_strand.",plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(coexpr ~ strands + sameTAD, data = all_ds_values, main = myTit,
        ylab = "Gene coexpression (PCC)")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# table(coexprDT$strands, coexprDT$sameTAD)


outFile <- file.path(outFold, paste0("all_datasets", "_multiline_density_coexpr_TAD_strand.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot_multiDens(list(
  divergent_sameTAD = all_ds_values$coexpr[all_ds_values$sameTAD == "sameTAD" & all_ds_values$strands == "diffStrand"],
  same_sameTAD = all_ds_values$coexpr[all_ds_values$sameTAD == "sameTAD" & all_ds_values$strands == "sameStrand"],
  divergent_diffTAD = all_ds_values$coexpr[all_ds_values$sameTAD == "diffTAD" & all_ds_values$strands == "diffStrand"],
  same_diffTAD = all_ds_values$coexpr[all_ds_values$sameTAD == "diffTAD" & all_ds_values$strands == "sameStrand"]
), my_xlab =  "Gene coexpression (PCC)", plotTit  = myTit, legPos = "topleft"
)
mtext(side = 3,text = paste0("(n = ", nDS, " datasets)"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


