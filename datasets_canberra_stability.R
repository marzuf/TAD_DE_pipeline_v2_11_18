startTime <- Sys.time()
cat(paste0("> START datasets_canberra_stability.R\n"))

# Rscript datasets_canberra_stability.R

source("analysis_utils.R")
options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

vdHeight <- 7
vdWidth <- 7

pvalSelect <- 0.1
withPvalSelect <- TRUE

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

outFold <- file.path("DATASET_CANBERRA_STABILITY")
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
    
    setNames(match(ref_tads, names(tad_pval)), ref_tads)

  } # end iterating datasets

  # rownames(all_ds_TAD_ranks_DT) <- all_ds
      
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


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

