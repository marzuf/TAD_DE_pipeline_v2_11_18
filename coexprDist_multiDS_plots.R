options(scipen=100)

library(foreach)
library(doMC)
library(tools)
library(Hmisc)

startTime <- Sys.time()

# Rscript coexprDist_multiDS_plots.R 
# Rscript coexprDist_multiDS_plots.R 10

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

outFold <- file.path("COEXPRDIST_MULTIDS_PLOTS", nTopDS_toPlot)
system(paste0("mkdir -p ", outFold))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) >= nTopDS_toPlot)


##########################################################################################
########################################################################################## RETRIEVE FCC
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
##########################################################################################
##########################################################################################

cat("... start building DT for the datasets\n")
all_coexpr_vect <- foreach(ds = all_ds) %dopar% {
  # all_coexpr_vect <- foreach(ds = all_ds, .combine='cbind') %dopar% {
    dist_vect_file <- file.path(mainDir,ds,
                         "distVect.Rdata")  
  stopifnot(file.exists(dist_vect_file))
  dist_vect <- eval(parse(text = load(dist_vect_file)))
  
  
  sameTADvect_file <- file.path(mainDir,ds,
                              "smooth_vals_sameTAD_distVect.Rdata")  
  stopifnot(file.exists(sameTADvect_file))
  sameTAD_vect <- eval(parse(text = load(sameTADvect_file)))
  
  
  diffTADvect_file <- file.path(mainDir,ds,
                                "smooth_vals_diffTAD_distVect.Rdata")  
  stopifnot(file.exists(diffTADvect_file))
  diffTAD_vect <- eval(parse(text = load(diffTADvect_file)))
  
  # data.frame(dist_vect = dist_vect,
    list(dist_vect = dist_vect,
              sameTAD_vect = sameTAD_vect,
              diffTAD_vect = diffTAD_vect,
             diff_vect = sameTAD_vect-diffTAD_vect)
}
names(all_coexpr_vect) <- all_ds
# colnames(all_coexpr_vect) <- paste0(
#                                     rep(all_ds, each=4),
#                                     "_",
#                                     colnames(all_coexpr_vect))

    # outFile <- file.path(outFold, "all_coexpr_vect.Rdata")
    # save(all_coexpr_vect, file = outFile)
    # cat(paste0("... written: ", outFile, "\n"))
    # stop("--ok--\n")
    # all_coexpr_vect <- eval(parse(text = load("COEXPRDIST_MULTIDS_PLOT_RESCALED/all_coexpr_vect.Rdata")))
    # str(all_coexpr_vect)
# all_coexpr_vect[1:5,1:5]
# check_dist <- all_coexpr_vect[,grepl("dist_vect", colnames(all_coexpr_vect))]
# stopifnot(nrow(unique(t(check_dist))) == 1)
stopifnot(length(unique(lapply(all_coexpr_vect, function(x) x[["dist_vect"]]))) == 1)


#==================================================================== PREPARE VARIABLES FOR THE PLOTS

cat("... start preparing values to plot\n")

distVect <- unlist(unique(lapply(all_coexpr_vect, function(x) x[["dist_vect"]])))

sameTAD_list <- lapply(all_coexpr_vect, function(x) x[["sameTAD_vect"]])
diffTAD_list <- lapply(all_coexpr_vect, function(x) x[["diffTAD_vect"]])

stopifnot(unique(unlist(lapply(diffTAD_list, length))) == unique(unlist(lapply(sameTAD_list, length))))
stopifnot(unique(unlist(lapply(diffTAD_list, length))) == length(distVect))

stopifnot(length(sameTAD_list) == length(diffTAD_list))

myxlab <- "Distance between gene pairs (bp)"
myylab <- "Pairwise coexpression (PCC qqnorm data)"

myTit <- paste0("Pairwise coexpr. ~ dist.")

mySub <- paste0("all datasets (n=", length(sameTAD_list), ")")

#==================================================================== START PLOTTING

### CURVES FOR DIFFERENT/SAME TADs

for(tadType in c("diff", "same")) {
  
  curr_list <- eval(parse(text = paste0(tadType, "TAD_list")))

  outFile <- file.path(outFold, paste0(tadType, "TAD_curves_allDS.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(NULL,
       xlim = range(distVect),
       # ylim = diffTAD_range,
       ylim = range(unlist(curr_list), na.rm=T),
       xlab = myxlab,
       ylab = myylab,
       main = paste0(myTit, " - ", tadType, "TAD")
  )
  mtext(side = 3, text = mySub)   
  
  foo <- lapply(1:length(curr_list), function(i) {
    dsname <- names(curr_list)[i]
    x <- unlist(curr_list[i])
    lineCol <- ifelse(dsname %in% topDS, topCol,
                      ifelse(dsname %in% botDS, botCol, otherCol))
    lines(x=distVect, y=x,
          col = lineCol)
  })
  
  foo <- lapply(1:length(curr_list), function(i) {
    dsname <- names(curr_list)[i]
    x <- unlist(curr_list[i])
    lineCol <- ifelse(dsname %in% topDS, topCol,
                      ifelse(dsname %in% botDS, botCol, otherCol))
    if(dsname %in% c(topDS, botDS)){
      lines(x=distVect, y=x,
            col = lineCol)
    }
  })
  
  legend("topright", lty=1, legend = c(paste0("top ", nTopDS_toPlot, " FCC AUC ratio DS"), paste0("last ", nTopDS_toPlot, " FCC AUC ratio DS")),
         col=c(topCol, botCol), bty="n")
  foo <- dev.off()
  
}

myylab <- "Diff. pairwise coexpression (PCC qqnorm data)"

curr_list <- lapply(all_coexpr_vect, function(x) x[["diff_vect"]])

outFile <- file.path(outFold, paste0("sameTAD_diffTAD_diff_curves_allDS.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     xlim = range(distVect),
     # ylim = diffTAD_range,
     ylim = range(unlist(curr_list), na.rm=T),
     xlab = myxlab,
     ylab = myylab,
     main = paste0(myTit, " - diff. sameTAD-diffTAD")
)
mtext(side = 3, text = mySub)   

foo <- lapply(1:length(curr_list), function(i) {
  dsname <- names(curr_list)[i]
  x <- unlist(curr_list[i])
  lineCol <- ifelse(dsname %in% topDS, topCol,
                    ifelse(dsname %in% botDS, botCol, otherCol))
  lines(x=distVect, y=x,
        col = lineCol)
})

foo <- lapply(1:length(curr_list), function(i) {
  dsname <- names(curr_list)[i]
  x <- unlist(curr_list[i])
  lineCol <- ifelse(dsname %in% topDS, topCol,
                    ifelse(dsname %in% botDS, botCol, otherCol))
  if(dsname %in% c(topDS, botDS)){
    lines(x=distVect, y=x,
          col = lineCol)
  }
})

abline(h=0, lty=2, lwd=0.5, col="red")

legend("topright", lty=1, legend = c(paste0("top ", nTopDS_toPlot, " FCC AUC ratio DS"), paste0("last ", nTopDS_toPlot, " FCC AUC ratio DS")),
       col=c(topCol, botCol), bty="n")
foo <- dev.off()

  
################################################################################
################################################################################
################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


