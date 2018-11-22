options(scipen=100)

library(foreach)
library(doMC)
library(tools)
library(Hmisc)

startTime <- Sys.time()

# Rscript fcc_multiDS_plots.R 
# Rscript fcc_multiDS_plots.R 10

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 30))

nTopDS_toPlot <- 5 # [default value]
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1){
  if(!is.na(as.numeric(args[1])))
    nTopDS_toPlot <- as.numeric(args[1])
}
cat(paste0("... nTopDS_toPlot =\t", nTopDS_toPlot, "\n"))

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)

my_xlab <- "Distance between gene pairs (bp)"
my_ylab <- "Gene pairs expression (qqnorm) PCC"
# densityPolygon <- 90
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 350, 8)

topCol <-  "darkorange2"
botCol <-  "dodgerblue4"
otherCol <- "gainsboro"

mainDir <- file.path(setDir,
                     "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                     "AUC_COEXPRDIST_SORTNODUP/")

caller <- "TopDom"

outFold <- file.path("FCC_MULTIDS_PLOTS", nTopDS_toPlot)
system(paste0("mkdir -p ", outFold))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) >= nTopDS_toPlot)


##########################################################################################
########################################################################################## RETRIEVE FCC VECTOR
##########################################################################################


vectFile <- file.path("WAVE_PLOT_AUCFCC", "all_plotVect.Rdata")
stopifnot(file.exists(vectFile))

all_fcc_vect <- eval(parse(text = load(vectFile)))


##########################################################################################
########################################################################################## RETRIEVE FCC AUC RATIO 
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

topDS <- names(sort(all_ds_aucFCC, decreasing=TRUE)[1:nTopDS_toPlot])
botDS <- names(rev(sort(all_ds_aucFCC, decreasing=TRUE))[1:nTopDS_toPlot])

cat(paste0("... will draw lines for ", length(all_ds) , " datasets\n"))

##########################################################################################
########################################################################################## PREPARE FOR THE PLOT
##########################################################################################


# filter to take the smallest number of TADs
min_TADnbr <- min(unlist(lapply(all_fcc_vect, length)))
stopifnot(!is.na(min_TADnbr))
stopifnot(min_TADnbr > 0)

all_fcc_vect_filter <- lapply(all_fcc_vect, function(x) x[1:min_TADnbr])
stopifnot(length(unique(lapply(all_fcc_vect_filter, length))) == 1)
stopifnot(unique(lapply(all_fcc_vect_filter, length)) == min_TADnbr)

xVect <- 1:min_TADnbr

myxlab <- "TAD rank"
myylab <- "cumsum(abs(FCC score))"

myTit <- paste0("FCC score cumsum")

mySub <- paste0("all datasets (n=", length(all_fcc_vect), "; min # TADs=", min_TADnbr, ")" )

#==================================================================== START PLOTTING

curr_list <- all_fcc_vect_filter

outFile <- file.path(outFold, paste0("wave_plots_allDS.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     xlim = range(xVect),
     ylim = range(unlist(curr_list), na.rm=T),
     xlab = myxlab,
     ylab = myylab,
     main = paste0(myTit)
)
mtext(side = 3, text = mySub)   

foo <- lapply(1:length(curr_list), function(i) {
  dsname <- names(curr_list)[i]
  x <- unlist(curr_list[i])
  lineCol <- ifelse(dsname %in% topDS, topCol,
                    ifelse(dsname %in% botDS, botCol, otherCol))
  lines(x=xVect, y=x,
        col = lineCol)
})

foo <- lapply(1:length(curr_list), function(i) {
  dsname <- names(curr_list)[i]
  x <- unlist(curr_list[i])
  lineCol <- ifelse(dsname %in% topDS, topCol,
                    ifelse(dsname %in% botDS, botCol, otherCol))
  if(dsname %in% c(topDS, botDS)){
    lines(x=xVect, y=x,
          col = lineCol)
  }
})


legend("topleft", lty=1, legend = c(paste0("top ", nTopDS_toPlot, " FCC AUC ratio DS"), paste0("last ", nTopDS_toPlot, " FCC AUC ratio DS")),
       col=c(topCol, botCol), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
  
################################################################################
################################################################################
################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


