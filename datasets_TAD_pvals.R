startTime <- Sys.time()
cat(paste0("> START datasets_TAD_pvals.R\n"))

# Rscript datasets_TAD_pvals.R

source("analysis_utils.R")
options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight*1.2

plotCex <- 1.2

vdHeight <- 7
vdWidth <- 7

signifThresh <- 0.05


source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18")


# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
# file with coordinates of all regions
TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
tad_DT <- read.delim(TADpos_file, stringsAsFactors = F, header=F, col.names=c("chromo", "region", "start", "end"))
stopifnot(is.numeric(tad_DT$start))
tad_DT <- tad_DT[grep("_TAD", tad_DT$region),]
tad_DT <- tad_DT[order(tad_DT$chromo, tad_DT$start, tad_DT$end),]
head(tad_DT)  
ref_tads <- tad_DT$region

outFold <- file.path("DATASETS_TAD_PVALS")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "canberra_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) > 0)

#all_ds <- all_ds[all_ds %in% run1_DS]

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)


##########################################################################################
##########################################################################################


curr_ds="TCGAcoad_msi_mss"

topDS <- "TCGAcoad_msi_mss"

topDS <- all_ds[1:3]
topDS <- all_ds

topDS <- curr_ds

stopifnot(length(topDS) == 1)

if(buildTable){
  all_ds_TAD_pvals_DT <- foreach(curr_ds = all_ds, .combine='rbind') %do% {
    # all_ds_DT <- foreach(curr_ds = topDS, .combine='rbind') %dopar% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    printAndLog(txt, logFile)

    ### RETRIEVE TAD GENES AND PVALUES
    cat("... retrieve top TADs genes\n")
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(file.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    
    xx <- tad_pval[ref_tads]
    stopifnot(length(xx) == length(ref_tads))
    stopifnot(names(xx)[!is.na(xx)] == ref_tads[!is.na(xx)])
    names(xx) <- ref_tads
    xx

  } # end iterating datasets

  rownames(all_ds_TAD_pvals_DT) <- all_ds
      
  outFile <-     file.path(outFold, "all_ds_TAD_pvals_DT.Rdata")
  save(all_ds_TAD_pvals_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {# end if buildtable
    
  outFile <-     file.path(outFold, "all_ds_TAD_pvals_DT.Rdata")
  all_ds_TAD_pvals_DT <- eval(parse(text = load(outFile)))
  
}

# load("DATASETS_TAD_PVALS/all_ds_TAD_pvals_DT.Rdata")
all_ds_TAD_pvals_DT[1:5,1:5]

# for each of the TADs -> in how many dataset it is signif
nSignif_byTADs <- apply(all_ds_TAD_pvals_DT, 2, function(x) sum(na.omit(x) <= signifThresh ))
range(nSignif_byTADs)

allDS_nSignif_byTADs <- nSignif_byTADs

myTit <- paste0("# datasets in which TAD signif.")
myxlab <- "# datasets in which signif."
myylab <- "density"
mySub <- paste0("(signif. threshold = ", signifThresh, ")")
top5 <- names(sort(nSignif_byTADs, decreasing=T))[1:5]

legTxt <-  paste0("# datasets = ", nrow(all_ds_TAD_pvals_DT),"\n", "# TADs = ", ncol(all_ds_TAD_pvals_DT), "\n", "top5:\n", paste0(top5, collapse=",\n"))


outFile <- file.path(outFold, paste0("nbr_all_ds_in_which_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(nSignif_byTADs),
     main=myTit,
     xlab=myxlab,
     ylab=myylab,
     cex.axis=plotCex, cex.lab=plotCex
     )
mtext(side=3, text = mySub)
legend( "topright", bty="n",legend=legTxt, lty=-1)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##################
######### TCGA data only
##################

tcga_ds_TAD_pvals_DT <- all_ds_TAD_pvals_DT[grepl("^TCGA", rownames(all_ds_TAD_pvals_DT)),]

nSignif_byTADs <- apply(tcga_ds_TAD_pvals_DT, 2, function(x) sum(na.omit(x) <= signifThresh ))
range(nSignif_byTADs)

tcgaDS_nSignif_byTADs <- nSignif_byTADs

myTit <- paste0("# datasets in which TAD signif.")
myxlab <- "# datasets in which signif."
myylab <- "density"
mySub <- paste0("(TCGA data; signif. threshold = ", signifThresh, ")")
top5 <- names(sort(nSignif_byTADs, decreasing=T))[1:5]

legTxt <-  paste0("# datasets = ", nrow(tcga_ds_TAD_pvals_DT),"\n", "# TADs = ", ncol(tcga_ds_TAD_pvals_DT), "\n", "top5:\n", paste0(top5, collapse=",\n"))

outFile <- file.path(outFold, paste0("nbr_tcga_ds_in_which_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(nSignif_byTADs),
     main=myTit,
     xlab=myxlab,
     ylab=myylab,
     cex.axis=plotCex, cex.lab=plotCex
     )
mtext(side=3, text = mySub)
legend( "topright", bty="n",legend=legTxt, lty=-1)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


stopifnot(length(allDS_nSignif_byTADs) == length(tcgaDS_nSignif_byTADs))
outFile <- file.path(outFold, paste0("nbr_tcga_vs_all_ds_in_which_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
plot(x=allDS_nSignif_byTADs, 
     y=tcgaDS_nSignif_byTADs,
     cex.axis=plotCex, cex.lab=plotCex,
     pch=16, cex=0.7,
     main = "# signif. all vs. TCGA datasets",
     xlab=paste0("# datasets where signif. (all datasets, n=", nrow(all_ds_TAD_pvals_DT), ")"),
     ylab =paste0("# datasets where signif. (TCGA datasets, n=", nrow(tcga_ds_TAD_pvals_DT), ")")
     )
mtext(side=3, text=paste0("# TADs=", length(allDS_nSignif_byTADs)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

stopifnot(length(allDS_nSignif_byTADs) == length(tcgaDS_nSignif_byTADs))
outFile <- file.path(outFold, paste0("nbr_tcga_vs_all_ds_in_which_signif_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(x=allDS_nSignif_byTADs, 
     y=tcgaDS_nSignif_byTADs,
     cex.axis=plotCex, cex.lab=plotCex,
     pch=16, cex=0.7,
     main = "# signif. all vs. TCGA datasets",
     xlab=paste0("# datasets where signif. (all datasets, n=", nrow(all_ds_TAD_pvals_DT), ")"),
     ylab =paste0("# datasets where signif. (TCGA datasets, n=", nrow(tcga_ds_TAD_pvals_DT), ")")
)
mtext(side=3, text=paste0("# TADs=", length(allDS_nSignif_byTADs)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# 
# commonTADs_DT <- all_ds_TAD_pvals_DT[,apply(all_ds_TAD_pvals_DT, 2,function(x) !any(is.na(x)))]
# stopifnot(!any(is.na(commonTADs_DT)))
# 
# cat("dim(all_ds_TAD_pvals_DT)\t=\t",dim(all_ds_TAD_pvals_DT), "\n" )
# cat("dim(commonTADs_DT)\t=\t",dim(commonTADs_DT), "\n" )
# 
# source("canberra_clean.R")
# 
# x = matrix(c(2,4,1,3,0,3,4,1,2,0,2,4,3,0,1), byrow = T, ncol=5)
# canberra_stability(commonTADs_DT[1:3,1:5])
# tt <- commonTADs_DT[1:3,1:5]
# rank_tt <- t(apply(tt, 1, function(x) rank(x)))
# canberra_stability(x)
# canberra_stability(tt)
# canberra_stability(rank_tt)
# 
# rank_commonTADs_DT <- t(apply(commonTADs_DT, 1, function(x) rank(x)))
# canberra_stability(rank_commonTADs_DT[1:3,1:5])
# 
# # stability all datasets:
# # head(rank_commonTADs_DT)
# cat("... stability all datasets:\n")
# canberra_stability(rank_commonTADs_DT)
# 
# # stability TCGA
# cat("... stability TCGA datasets:\n")
# xx <- rank_commonTADs_DT[grepl("^TCGA", rownames(rank_commonTADs_DT)), ]
# dim(xx)
# canberra_stability(xx)
# 
# # stability no TCGA
# cat("... stability no TCGA datasets:\n")
# canberra_stability(rank_commonTADs_DT[! grepl("^TCGA", rownames(rank_commonTADs_DT)), ])
# 
# 
# 
# #### no pval select
# # dim(all_ds_TAD_pvals_DT)        =        84 3986 
# # dim(commonTADs_DT)      =        84 692 
# # [1] 0
# # [1] 0.6688182
# # [1] 0
# # [1] 0.416203
# # [1] 0
# # ... stability all datasets:
# #   [1] 0.9328141
# # ... stability TCGA datasets:
# #   [1] 0.8529246
# # ... stability no TCGA datasets:
# #   [1] 0.9513543
# 
# ### with pval select 0.05
# # load("DATASETS_TAD_PVALS/all_ds_TAD_pvals_DT.Rdata")
# all_ds_TAD_pvals_DT[1:5,1:5]
# 
# rankMax <- max(all_ds_TAD_pvals_DT, na.rm=T)
# 
# cat(paste0("... start drawing\n"))
# outFile <- file.path(outFold, paste0("all_TADs_cumsum_linePlot.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot(NULL, 
#      xlim=c(1,rankMax),
#      ylim=c(1,nrow(all_ds_TAD_pvals_DT)),
#      main=paste0("cumsum # datasets <= rank"),
#      ylab = paste0("# datasets"),
#      xlab = paste0("TAD rank")
#      )
# for(i in 1:ncol(all_ds_TAD_pvals_DT)) {
#   curr_tad_ranks <- all_ds_TAD_pvals_DT[,i]
#   rankVect <- 1:rankMax
#   cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
#   lines(x=rankVect, y= cumsum_ranks)
# }
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# cat(paste0("... start computing auc\n"))
# all_TADs_auc <- foreach(i = 1:ncol(all_ds_TAD_pvals_DT), .combine='c') %dopar% {
#   curr_tad_ranks <- all_ds_TAD_pvals_DT[,i]
#   rankVect <- 1:rankMax
#   cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
#   curr_auc <- auc(x = rankVect, y = cumsum_ranks)
#   curr_auc
# }
# outFile <- file.path(outFold, "all_TADs_auc.Rdata")
# save(all_TADs_auc, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# #### DO THE SAME FOR TCGA ONLY
# 
# all_ds_TAD_pvals_DT[1:5,1:5]
# 
# tcga_ds_TAD_ranks_DT <- all_ds_TAD_pvals_DT[grep("^TCGA", rownames(all_ds_TAD_pvals_DT)),]
# 
# rankMax <- max(tcga_ds_TAD_ranks_DT, na.rm=T)
# 
# cat(paste0("... start drawing\n"))
# outFile <- file.path(outFold, paste0("tcga_TADs_cumsum_linePlot.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot(NULL, 
#      xlim=c(1,rankMax),
#      ylim=c(1,nrow(tcga_ds_TAD_ranks_DT)),
#      main=paste0("cumsum # datasets <= rank"),
#      ylab = paste0("# datasets"),
#      xlab = paste0("TAD rank")
# )
# for(i in 1:ncol(tcga_ds_TAD_ranks_DT)) {
#   curr_tad_ranks <- tcga_ds_TAD_ranks_DT[,i]
#   rankVect <- 1:rankMax
#   cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
#   lines(x=rankVect, y= cumsum_ranks)
# }
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# cat(paste0("... start computing TCGA auc\n"))
# tcga_TADs_auc <- foreach(i = 1:ncol(tcga_ds_TAD_ranks_DT), .combine='c') %dopar% {
#   curr_tad_ranks <- tcga_ds_TAD_ranks_DT[,i]
#   rankVect <- 1:rankMax
#   cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
#   curr_auc <- auc(x = rankVect, y = cumsum_ranks)
#   curr_auc
# }
# outFile <- file.path(outFold, "tcga_TADs_auc.Rdata")
# save(tcga_TADs_auc, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

