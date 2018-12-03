# if promoter sharing as a source of coexpression,
# divergently transcriped genes should have higher coexpr.

setDir <- "/media/electron"
setDir <- ""

source("analysis_utils.R")

plotType <- "svg"
myHeight <- 7
myWidth <- 10

plotType <- "png"
myHeight <- 400
myWidth <- 600


outFold <- "COEXPR_TAD_STRAND"
system(paste0("mkdir -p ", outFold))

curr_dataset <- "TCGAcrc_msi_mss"
myTit <- paste0(curr_dataset, " - Coexpression: TAD and strand ")

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

outFile <- file.path(outFold, paste0(curr_dataset, "_boxplot_coexpr_TAD_strand.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(coexpr ~ strands + sameTAD, data = coexprDT, main = myTit,
        ylab = "Gene coexpression (PCC)")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

table(coexprDT$strands, coexprDT$sameTAD)


outFile <- file.path(outFold, paste0(curr_dataset, "_multiline_density_coexpr_TAD_strand.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot_multiDens(list(
  divergent_sameTAD = coexprDT$coexpr[coexprDT$sameTAD == "sameTAD" & coexprDT$strands == "diffStrand"],
  same_sameTAD = coexprDT$coexpr[coexprDT$sameTAD == "sameTAD" & coexprDT$strands == "sameStrand"],
  divergent_diffTAD = coexprDT$coexpr[coexprDT$sameTAD == "diffTAD" & coexprDT$strands == "diffStrand"],
  same_diffTAD = coexprDT$coexpr[coexprDT$sameTAD == "diffTAD" & coexprDT$strands == "sameStrand"]
), my_xlab =  "Gene coexpression (PCC)", plotTit  = myTit,  legPos = "topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
