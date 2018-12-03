library(foreach)
library(doMC)

# Rscript all_AUC_ranking.R
cat("> START: all_AUC_ranking.R\n")

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

startTime <- Sys.time()


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

outFold <- file.path(paste0("ALL_AUC_RANKING"))
system(paste0("mkdir -p ", outFold))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)


##########################################################################################
##########################################################################################
##########################################################################################

all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  stopifnot(file.exists(aucFCC_file))
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  aucFCC
}
names(all_ds_aucFCC) <- all_ds

all_ds_aucFCC <- sort(all_ds_aucFCC, decreasing=TRUE)

stopifnot(all_ds %in% names(all_ds_aucFCC))

##########################################################################################
##########################################################################################
##########################################################################################

all_ds_aucCoexpr <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  aucCoexprDist
}
names(all_ds_aucCoexpr) <- all_ds

all_ds_aucCoexpr <- sort(all_ds_aucCoexpr, decreasing=TRUE)

stopifnot(all_ds %in% names(all_ds_aucCoexpr))


##########################################################################################
##########################################################################################
##########################################################################################

all_ds_aucFC_file <- file.path("AUC_RATIO_CUMULVECT",
                                    "FC", "raw", "auc_values.Rdata")

all_ds_aucFC <- eval(parse(text = load(all_ds_aucFC_file)))
all_ds_aucFC <- unlist(all_ds_aucFC[["rescAUC"]])
all_ds_aucFC <- sort(all_ds_aucFC, decreasing=TRUE)

stopifnot(all_ds %in% names(all_ds_aucFC))

##########################################################################################
##########################################################################################
##########################################################################################

all_ds_aucVar_file <- file.path("AUC_RATIO_CUMULVECT", 
                                "var", "raw" ,
                                "auc_values.Rdata")

all_ds_aucVar <- eval(parse(text = load(all_ds_aucVar_file)))
all_ds_aucVar <- unlist(all_ds_aucVar[["rescAUC"]])
all_ds_aucVar <- sort(all_ds_aucVar, decreasing=TRUE)

stopifnot(all_ds %in% names(all_ds_aucVar))


####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################


intersectList <- list(
  all_vects =  c("all_ds_aucVar", "all_ds_aucFC", "all_ds_aucCoexpr", "all_ds_aucFCC"),
  fcc_coexpr =  c("all_ds_aucCoexpr", "all_ds_aucFCC")
)

# all_vects <- c("all_ds_aucVar", "all_ds_aucFC", "all_ds_aucCoexpr", "all_ds_aucFCC")


for(k in seq_along(intersectList)) {

  curr_inter <- intersectList[[k]]
    
  subTit <- paste0("vectors: ", paste0(gsub("all_ds_auc", "", curr_inter), collapse=", "))
  
  nAtIntersect <- foreach(i = c(1:length(all_ds)), .combine='c') %dopar% {
    all_tops <- lapply(curr_inter, function(x) {
      curr_vect <- eval(parse(text=x))
      names(curr_vect[1:i])
    })
    intersectDS <- Reduce(intersect, all_tops)
    length(intersectDS)
  }
  names(nAtIntersect) <- paste0(c(1:length(all_ds)))
  
  
  saveIntersect <- foreach(i = c(1:length(all_ds))) %dopar% {
    all_tops <- lapply(curr_inter, function(x) {
      curr_vect <- eval(parse(text=x))
      names(curr_vect[1:i])
    })
    intersectDS <- Reduce(intersect, all_tops)
    intersectDS
  }
  outFile <- file.path(outFold, paste0("saveIntersect_", names(intersectList)[k], ".Rdata"))
  save(saveIntersect, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  axLabs <- seq(from=1, to=length(all_ds), by=5)
  
  outFile <- file.path(outFold, paste0("intersectTopDs_by_nTop_", names(intersectList)[k], ".", plotType))
  
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(x = as.numeric(as.character(names(nAtIntersect))), 
       y = nAtIntersect,
       main = "Intersect top datasets",
       ylab = "# datasets at the intersect",
       xlab = "# top datasets",
       type="l",
       bty ="L",
       cex = 1.5,
       cex.lab=1.2, axes =F)
  axis(1, at = axLabs,labels = axLabs)
  axis(2, at = axLabs, labels = axLabs)
  mtext(side = 3, text = paste0(subTit))
  box(bty="L")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}



######################################################################################
######################################################################################
######################################################################################


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






