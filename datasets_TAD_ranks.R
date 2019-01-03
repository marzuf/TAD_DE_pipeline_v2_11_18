startTime <- Sys.time()
cat(paste0("> START datasets_TAD_ranks.R\n"))

# Rscript datasets_TAD_ranks.R

source("analysis_utils.R")
options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "png"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight*1.2

vdHeight <- 7
vdWidth <- 7

pvalSelect <- 0.1
withPvalSelect <- F

nTopToPlot = 5

computeAUC <- TRUE
cat("!!! compute AUC = ", as.character(computeAUC), "\n")

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

outFold <- file.path("DATASETS_TAD_RANKS")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "canberra_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) > 0)

#all_ds <- all_ds[all_ds %in% run1_DS]

txt <- paste0("... computeAUC:\t", as.character(computeAUC), "\n")
printAndLog(txt, logFile)

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
  all_ds_TAD_ranks_DT <- foreach(curr_ds = all_ds, .combine='rbind') %do% {
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
    tad_pval <- sort(p.adjust(tad_pval, method = "BH"))
    if(withPvalSelect){
      tad_pval <- tad_pval[tad_pval <= pvalSelect]
      if(length(tad_pval) < 1000) return(NULL)
    }
    tad_ranks <- rank(tad_pval, ties="min")
    # > tt
    # chr1_TAD150  chr1_TAD227 chr11_TAD143 
    # 0.005255629  0.005255629  0.005255629 
    # > ref <- c( "chr11_TAD143","chr1_TAD150", "chr1_TAD227")
    # > match(ref, names(tt))
    # [1] 3 1 2
    # > ref <- c(  "chr1_TAD227","chr11_TAD143","chr1_TAD150")
    # > match(ref, names(tt))
    # [1] 2 3 1
    # > ref <- c(  "chr1_TAD227","chr11_TAD143","chr1_TAD150", "chr1_TAD6")
    # > match(ref, names(tt))
    # [1]  2  3  1 NA
    
    # Canberra fct: row = dataset (position list), column = position
    
    # setNames(match(ref_tads, names(tad_pval)), ref_tads)
    
    xx <- tad_ranks[ref_tads]
    stopifnot(length(xx) == length(ref_tads))
    stopifnot(names(xx)[!is.na(xx)] == ref_tads[!is.na(xx)])
    names(xx) <- ref_tads
    xx

  } # end iterating datasets

  rownames(all_ds_TAD_ranks_DT) <- all_ds
      
  outFile <-     file.path(outFold, "all_ds_TAD_ranks_DT.Rdata")
  save(all_ds_TAD_ranks_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {# end if buildtable
    
  outFile <-     file.path(outFold, "all_ds_TAD_ranks_DT.Rdata")
  all_ds_TAD_ranks_DT <- eval(parse(text = load(outFile)))
  
}


commonTADs_DT <- all_ds_TAD_ranks_DT[,apply(all_ds_TAD_ranks_DT, 2,function(x) !any(is.na(x)))]
stopifnot(!any(is.na(commonTADs_DT)))

cat("dim(all_ds_TAD_ranks_DT)\t=\t",dim(all_ds_TAD_ranks_DT), "\n" )
cat("dim(commonTADs_DT)\t=\t",dim(commonTADs_DT), "\n" )

source("canberra_clean.R")

x = matrix(c(2,4,1,3,0,3,4,1,2,0,2,4,3,0,1), byrow = T, ncol=5)
canberra_stability(commonTADs_DT[1:3,1:5])
tt <- commonTADs_DT[1:3,1:5]
rank_tt <- t(apply(tt, 1, function(x) rank(x)))
canberra_stability(x)
canberra_stability(tt)
canberra_stability(rank_tt)

rank_commonTADs_DT <- t(apply(commonTADs_DT, 1, function(x) rank(x)))
canberra_stability(rank_commonTADs_DT[1:3,1:5])

# stability all datasets:
# head(rank_commonTADs_DT)
cat("... stability all datasets:\n")
canberra_stability(rank_commonTADs_DT)

# stability TCGA
cat("... stability TCGA datasets:\n")
xx <- rank_commonTADs_DT[grepl("^TCGA", rownames(rank_commonTADs_DT)), ]
dim(xx)
canberra_stability(xx)

# stability no TCGA
cat("... stability no TCGA datasets:\n")
canberra_stability(rank_commonTADs_DT[! grepl("^TCGA", rownames(rank_commonTADs_DT)), ])



#### no pval select
# dim(all_ds_TAD_ranks_DT)        =        84 3986 
# dim(commonTADs_DT)      =        84 692 
# [1] 0
# [1] 0.6688182
# [1] 0
# [1] 0.416203
# [1] 0
# ... stability all datasets:
#   [1] 0.9328141
# ... stability TCGA datasets:
#   [1] 0.8529246
# ... stability no TCGA datasets:
#   [1] 0.9513543

### with pval select 0.05
# load("DATASETS_TAD_RANKS/all_ds_TAD_ranks_DT.Rdata")
all_ds_TAD_ranks_DT[1:5,1:5]

rankMax <- max(all_ds_TAD_ranks_DT, na.rm=T)

#### DRAW THE CURVES

cat(paste0("... start drawing\n"))
outFile <- file.path(outFold, paste0("all_TADs_cumsum_linePlot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL, 
     xlim=c(1,rankMax),
     ylim=c(1,nrow(all_ds_TAD_ranks_DT)),
     main=paste0("cumsum # datasets <= rank"),
     ylab = paste0("# datasets"),
     xlab = paste0("TAD rank")
     )
for(i in 1:ncol(all_ds_TAD_ranks_DT)) {
  curr_tad_ranks <- all_ds_TAD_ranks_DT[,i]
  rankVect <- 1:rankMax
  cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
  lines(x=rankVect, y= cumsum_ranks)
}
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

if(computeAUC) {
  cat(paste0("... start computing auc\n"))
  all_TADs_auc <- foreach(i = 1:ncol(all_ds_TAD_ranks_DT), .combine='c') %dopar% {
    curr_tad_ranks <- all_ds_TAD_ranks_DT[,i]
    rankVect <- 1:rankMax
    cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
    curr_auc <- auc(x = rankVect, y = cumsum_ranks)
    curr_auc
  }
  names(all_TADs_auc) <- colnames(all_ds_TAD_ranks_DT)
  outFile <- file.path(outFold, "all_TADs_auc.Rdata")
  save(all_TADs_auc, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_TADs_auc.Rdata")
  all_TADs_auc <- eval(parse(text = load(outFile)))
}

############ DRAW ONLY TOP TADS

stopifnot(nTopToPlot <= length(all_TADs_auc))
topTADs <- names(sort(all_TADs_auc, decreasing = TRUE)[1:nTopToPlot])
stopifnot(topTADs %in% colnames(all_ds_TAD_ranks_DT) )
top_all_ds_TAD_ranks_DT <- all_ds_TAD_ranks_DT[, topTADs]

cat(paste0("... start drawing - only top AUC TADs \n"))
outFile <- file.path(outFold, paste0("all_TADs_cumsum_linePlot_AUCtop", nTopToPlot, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL, 
     xlim=c(1,rankMax),
     ylim=c(1,nrow(top_all_ds_TAD_ranks_DT)),
     main=paste0("cumsum # datasets <= rank"),
     ylab = paste0("# datasets"),
     xlab = paste0("TAD rank")
)
for(i in 1:ncol(top_all_ds_TAD_ranks_DT)) {
  curr_tad_ranks <- top_all_ds_TAD_ranks_DT[,i]
  rankVect <- 1:rankMax
  cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
  lines(x=rankVect, y= cumsum_ranks, col=i)
}

legend("topleft", lty=-1,bty="n", cex=0.7, legend = topTADs, col=1:length(topTADs)  )

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#### DO THE SAME FOR TCGA ONLY

all_ds_TAD_ranks_DT[1:5,1:5]

tcga_ds_TAD_ranks_DT <- all_ds_TAD_ranks_DT[grep("^TCGA", rownames(all_ds_TAD_ranks_DT)),]

rankMax <- max(tcga_ds_TAD_ranks_DT, na.rm=T)

cat(paste0("... start drawing\n"))
outFile <- file.path(outFold, paste0("tcga_TADs_cumsum_linePlot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL, 
     xlim=c(1,rankMax),
     ylim=c(1,nrow(tcga_ds_TAD_ranks_DT)),
     main=paste0("cumsum # datasets <= rank"),
     ylab = paste0("# datasets"),
     xlab = paste0("TAD rank")
)
for(i in 1:ncol(tcga_ds_TAD_ranks_DT)) {
  curr_tad_ranks <- tcga_ds_TAD_ranks_DT[,i]
  rankVect <- 1:rankMax
  cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
  lines(x=rankVect, y= cumsum_ranks)
}
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

if(computeAUC) {
  cat(paste0("... start computing TCGA auc\n"))
  tcga_TADs_auc <- foreach(i = 1:ncol(tcga_ds_TAD_ranks_DT), .combine='c') %dopar% {
    curr_tad_ranks <- tcga_ds_TAD_ranks_DT[,i]
    rankVect <- 1:rankMax
    cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
    curr_auc <- auc(x = rankVect, y = cumsum_ranks)
    curr_auc
  }
  outFile <- file.path(outFold, "tcga_TADs_auc.Rdata")
  save(tcga_TADs_auc, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "tcga_TADs_auc.Rdata")
  tcga_TADs_auc <- eval(parse(text = load(outFile)))
}


############ DRAW ONLY TOP TADS - TCGA ONLY

stopifnot(nTopToPlot <= length(tcga_TADs_auc))
tcga_ topTADs <- names(sort(tcga_TADs_auc, decreasing = TRUE)[1:nTopToPlot])
stopifnot(tcga_ topTADs %in% colnames(tcga_ds_TAD_ranks_DT) )
top_tcga_ds_TAD_ranks_DT <- tcga_ds_TAD_ranks_DT[, tcga_ topTADs]

cat(paste0("... start drawing - only top AUC TADs \n"))
outFile <- file.path(outFold, paste0("tcga_TADs_cumsum_linePlot_AUCtop", nTopToPlot, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL, 
     xlim=c(1,rankMax),
     ylim=c(1,nrow(top_tcga_ds_TAD_ranks_DT)),
     main=paste0("cumsum # datasets <= rank"),
     ylab = paste0("# datasets"),
     xlab = paste0("TAD rank")
)
for(i in 1:ncol(top_tcga_ds_TAD_ranks_DT)) {
  curr_tad_ranks <- top_tcga_ds_TAD_ranks_DT[,i]
  rankVect <- 1:rankMax
  cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_tad_ranks) <= x  ))
  lines(x=rankVect, y= cumsum_ranks, col=i)
}

legend("topleft", lty=-1,bty="n", cex=0.7, legend = tcga_ topTADs, col=1:length(tcga_ topTADs)  )

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

