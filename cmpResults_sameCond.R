startTime <- Sys.time()

library(foreach)
library(doMC)

source("analysis_utils.R")

cat("> START: cmpResults_sameCond.R\n")
# Rscript cmpResults_sameCond.R

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

pipFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER")

outFold <- "CMPRESULTS_SAMECOND"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmpResults_sameCond_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

pvalThresh <- 0.05

gene2tadFile <- file.path(setDir, 
                       "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 

gene2tad_DT <- read.delim(gene2tadFile, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)
stopifnot(!duplicated(gene2tad_DT$entrezID))
gene2tad_vect <- setNames(gene2tad_DT$region, gene2tad_DT$entrezID)

# /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/TCGAbrca_lum_bas/11_runEmpPvalCombined/emp_pval_combined.Rdata

all_similarCond <- list(
  c("GSE65540_before_after","GSE66306_before_after"),
  c("GSE58135_ERpos_tripleNeg", "TCGAbrca_lum_bas"),
  c("GSE81046_noninf_list", "GSE73765_noninf_list"),
  c("GSE73765_noninf_salm", "GSE81046_noninf_salm"),
  c("GSE73765_salm_list", "GSE81046_salm_list")
)


i=1
for(i in seq_along(all_similarCond)) {
  
  all_ds <- all_similarCond[[i]]
  
  ds = "TCGAbrca_lum_bas"
  all_ds_TADs <- foreach(ds = all_ds) %do% {
    
    
    ### PREPARE TADs
    tadpvalsFile <- file.path(pipFold, ds, "11_runEmpPvalCombined", "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalsFile))
    
    tadpvals <- eval(parse(text = load(tadpvalsFile)))
    tadpvals <- p.adjust(tadpvals, method="BH")
    tadpvals <- sort(tadpvals, decreasing=FALSE)
    signifTADs <- names(tadpvals[tadpvals <= pvalThresh])
    
    tadranks <- rank(tadpvals, ties = "min")
    
    ### PREPARE GENES
    geneListFile <- file.path(pipFold, ds, "0_prepGeneData", "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    geneList <- eval(parse(text = load(geneListFile)))
    
    gene2tads <- gene2tad_vect[names(gene2tad_vect) %in% geneList]    
    
    stopifnot(length(gene2tads) == length(geneList) )
    
    stopifnot(gene2tads %in% names(tadranks))
    
    genetadrank <- setNames(tadranks[gene2tads], names(gene2tads))

    txt <- paste0("... ", ds, " # TADs\t=\t", length(tadranks), "\n")
    printAndLog(txt, logFile)
    txt <- paste0("... ", ds, " # signif. TADs\t=\t", length(signifTADs), "\n")
    printAndLog(txt, logFile)
    txt <- paste0("... ", ds, " # genes\t=\t", length(geneList), "\n")
    printAndLog(txt, logFile)
    
    list(signifTADs=signifTADs, tadranks=tadranks, genetadrank = genetadrank)
  }
  
  nSignifTADs <- lapply(all_ds_TADs, function(x) length(x[["signifTADs"]]))
  nIntersectTADs <- length(Reduce(intersect, lapply(all_ds_TADs, function(x) x[["signifTADs"]])))
  
  nGenes <- lapply(all_ds_TADs, function(x) length(x[["genetadrank"]]))
  
  commonGenes <- Reduce(intersect, lapply(all_ds_TADs, function(x) names(x[["genetadrank"]])))
  stopifnot(commonGenes %in% gene2tad_DT$entrezID)
  
  nIntersectGenes <- length(commonGenes)
  
  all_gene2tad_ranks <- lapply(all_ds_TADs, function(x) x[["genetadrank"]][names(x[["genetadrank"]]) %in% commonGenes])

  stopifnot(length(Reduce(setdiff, lapply(all_gene2tad_ranks, names))) == 0)
  stopifnot(length(unique(lapply(all_gene2tad_ranks, names))) == 1)
  
  all_gene2tad_DT <- setNames(data.frame(all_gene2tad_ranks),all_ds)
  
  myx <- all_gene2tad_DT[,1]
  myy <- all_gene2tad_DT[,2]
  
  exprDS <- paste0(all_ds, collapse="_")
  
  myTit <- paste0(exprDS, ": gene TAD ranking comparison")
  myxlab <- paste0("gene TAD rank ", colnames(all_gene2tad_DT)[1], " (", nIntersectGenes, "/", nGenes[[1]], ")")
  myylab <- paste0("gene TAD rank ", colnames(all_gene2tad_DT)[2], " (", nIntersectGenes, "/", nGenes[[2]], ")")
  mySub <- paste0(all_ds[1], " vs. ", all_ds[2])
  
  outFile <- file.path(outFold, paste0(exprDS, "_gene_tad_rank_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  densplot(x=myx,
           y=myy,
           xlab=myxlab,
           ylab=myylab,
           pch = 16, cex = 0.7,
           main = myTit
  )
  mtext(side=3, text = mySub)
  add_curv_fit(x = myx,
               y= myy, withR2 = FALSE, lty=2, col="darkgray")
  
  addCorr(x=myx,
          y = myy,
          bty="n",
          legPos="bottomright")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
} 





######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



