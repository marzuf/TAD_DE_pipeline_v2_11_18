startTime <- Sys.time()

library(foreach)
library(doMC)

source("analysis_utils.R")

cat("> START: cmp_TAD_pvals.R\n")
# Rscript cmp_TAD_pvals.R

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

pipFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER")

outFold <- "CMP_TAD_PVALS"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_TAD_pvals_logFile.txt")
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
    
    txt <- paste0("... ", ds, " # TADs\t=\t", length(tadpvals), "\n")
    printAndLog(txt, logFile)
    txt <- paste0("... ", ds, " # signif. TADs\t=\t", length(signifTADs), "\n")
    printAndLog(txt, logFile)
    
    list(signifTADs=signifTADs, tadpvals=tadpvals)
  }
  
  ################################################# ALL TADs 
  intersectTADs <- Reduce(intersect, lapply(all_ds_TADs, function(x) names(x[["tadpvals"]])))
  
  length(all_ds_TADs[[1]][["tadpvals"]])
  length(all_ds_TADs[[2]][["tadpvals"]])
  length(intersectTADs)
  
  all_tad_pvals <- lapply(all_ds_TADs, function(x) x[["tadpvals"]][names(x[["tadpvals"]]) %in% intersectTADs])
  all_tad_pvals <- lapply(all_tad_pvals, function(x) x[intersectTADs])
  
  length(all_tad_pvals[[1]])
  length(all_tad_pvals[[2]])
  
  nTADs <- lapply(all_ds_TADs, function(x) length(x[["tadpvals"]]))
  nIntersectTADs <- length(intersectTADs)
  
  stopifnot(length(Reduce(setdiff, lapply(all_tad_pvals, names))) == 0)
  stopifnot(length(unique(lapply(all_tad_pvals, names))) == 1)
  
  all_tad_DT <- setNames(data.frame(all_tad_pvals),all_ds)
  
  myx <- -log10(all_tad_DT[,1])
  myy <- -log10(all_tad_DT[,2])
  
  exprDS <- paste0(all_ds, collapse="_")
  
  myTit <- paste0(exprDS, ": TAD pval comparison")
  myxlab <- paste0("-log10 TAD pval ", colnames(all_tad_DT)[1], " (", nIntersectTADs, "/", nTADs[[1]], ")")
  myylab <- paste0("-log10 TAD pval ", colnames(all_tad_DT)[2], " (", nIntersectTADs, "/", nTADs[[2]], ")")
  mySub <- paste0(all_ds[1], " vs. ", all_ds[2])
  
  if(length(myx) > 1) {
    outFile <- file.path(outFold, paste0(exprDS, "_tad_pval_density.", plotType))
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
  
  ################################################# SIGNIF TADs ONLY
  all_signifTAD_DT <- all_tad_DT[ all_tad_DT[,1] <= pvalThresh & all_tad_DT[,2] <= pvalThresh,  ]
  
  myx <- -log10(all_signifTAD_DT[,1])
  myy <- -log10(all_signifTAD_DT[,2])
  
  nintersectSignifTADs <- nrow(all_signifTAD_DT)
  
  exprDS <- paste0(all_ds, collapse="_")
  
  myTit <- paste0(exprDS, ": TAD pval comparison (signif. TADs only)")
  myxlab <- paste0("-log10 TAD pval ", colnames(all_signifTAD_DT)[1], " (", nintersectSignifTADs, "/", nTADs[[1]], ")")
  myylab <- paste0("-log10 TAD pval ", colnames(all_signifTAD_DT)[2], " (", nintersectSignifTADs, "/", nTADs[[2]], ")")
  mySub <- paste0(all_ds[1], " vs. ", all_ds[2])
  
  if(length(myx) > 1) {
    outFile <- file.path(outFold, paste0(exprDS, "_signifTAD_pval_density.", plotType))
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
  
} 




######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



