
# Rscript check_FPKM_geneLength_RSEMdatasets.R

add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}


source("analysis_utils.R")


plotType <- "png"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight


curr_ds <- "TCGAstad_msi_gs"
curr_ds <- "TCGAbrca_lum_bas"
curr_ds <- "TCGAucec_msi_cnl"

setDir <- "/media/electron"
setDir <- ""

outFold <- "CHECK_FPKM_GENE_LENGTH_RSEM"
system(paste0("mkdir -p ", outFold))

fpkm_v0 <- eval(parse(text = load(file.path(setDir,
                                            "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER_TCGA_RUN1",
                  curr_ds, "0_prepGeneData/rna_fpkmDT.Rdata"))))
fpkm_v0_geneMeanExpr <- rowMeans(fpkm_v0)

fpkm_v1 <- eval(parse(text = load(file.path(setDir,
                                            "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER",
                                curr_ds, "0_prepGeneData/rna_fpkmDT.Rdata"))))
fpkm_v1_geneMeanExpr <- rowMeans(fpkm_v1)

# HARD CODED
gffDT_file <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gffDT <- read.delim(gffDT_file, header=T, stringsAsFactors = F)
gffDT$entrezID <- as.character(gffDT$entrezID)

geneLength <- setNames(gffDT$end-gffDT$start, gffDT$entrezID)

# all(names(fpkm_v1_geneMeanExpr) %in% gffDT$entrezID)
# all(names(fpkm_v0_geneMeanExpr) %in% gffDT$entrezID)


fpkm_v0_geneMeanExpr <- fpkm_v0_geneMeanExpr[names(fpkm_v0_geneMeanExpr) %in% names(geneLength)]
head(fpkm_v0_geneMeanExpr)

# plot(fpkm_v0_geneMeanExpr ~ geneLength[names(fpkm_v0_geneMeanExpr)])

myx = log10(geneLength[names(fpkm_v0_geneMeanExpr)])
myy =log10(fpkm_v0_geneMeanExpr)
myx2 <- myx[!is.na(myx) & !is.na(myy) & !is.infinite(myx) & ! is.infinite(myy)]
myy2 <- myy[!is.na(myx) & !is.na(myy) & !is.infinite(myx) & ! is.infinite(myy)]
myx <- myx2
myy <- myy2

outFile <- file.path(outFold, paste0(curr_ds, "_wrong_FPKM.", plotType))
do.call(plotType, list(outFile, height= myHeight, width = myWidth))
plot( myy~ myx,
     xlab="log10 gene length",
     ylab = "log10 mean expr",
     main = paste0(curr_ds, " - wrong FPKM (expected counts)"),
     pch = 16, cex = 0.7, cex.axis = 1.2, cex.lab = 1.2
     )
add_curv_fit(x=myx, 
             y=myy
               )
addCorr(x=myx,
        y=myy
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

fpkm_v1_geneMeanExpr <- fpkm_v1_geneMeanExpr[names(fpkm_v1_geneMeanExpr) %in% names(geneLength)]
head(fpkm_v1_geneMeanExpr)

myx = log10(geneLength[names(fpkm_v1_geneMeanExpr)])
myy = log10(fpkm_v1_geneMeanExpr)
myx2 <- myx[!is.na(myx) & !is.na(myy) & !is.infinite(myx) & ! is.infinite(myy)]
myy2 <- myy[!is.na(myx) & !is.na(myy) & !is.infinite(myx) & ! is.infinite(myy)]
myx <- myx2
myy <- myy2


# plot(fpkm_v1_geneMeanExpr ~ geneLength[names(fpkm_v1_geneMeanExpr)])
outFile <- file.path(outFold, paste0(curr_ds, "_correct_FPKM.", plotType))
do.call(plotType, list(outFile, height= myHeight, width = myWidth))
plot(myy ~ myx,
     xlab="log10 gene length",
     ylab = "log10 mean expr",
     main = paste0(curr_ds, " - correct FPKM (scaled estimates)"),
     pch = 16, cex = 0.7, cex.axis = 1.2, cex.lab = 1.2
     )
add_curv_fit(x=myx, 
             y=myy
)
addCorr(x=myx,
        y=myy
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))