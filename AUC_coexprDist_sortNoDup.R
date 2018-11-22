startTime <- Sys.time()
cat(paste0("> Rscript AUC_coexprDist_sortNoDup.R\n"))

options(scipen=100)

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggstatsplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

### HARD CODED
caller <- "TopDom"
corMethod <- "pearson"
buildTable <- FALSE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"

# nbr of points for loess fit to take the AUC
nbrLoessPoints <- 1000

# UPDATE 30.06.2018:
# -> check that always $gene1 < $gene2 before left_join !!!

### RETRIEVE FROM COMMAND LINE
# Rscript coexpr_dist_v3.R

# top-ranking:
# Rscript AUC_coexprDist.R TCGAcrc_msi_mss

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 1)
  curr_dataset <- args[1]
} else{
  txt <- paste0("> Default parameters:\n")
  curr_dataset <- "TCGAcrc_msi_mss"
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "AUC_COEXPRDIST_SORTNODUP",  paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("AUC_coexprDist_logFile.txt"))  
system(paste0("rm -f ", logFile))

printAndLog(txt, logFile)

txt <- paste0("... curr_dataset = ",  curr_dataset, "\n")
printAndLog(txt, logFile)

txt <- paste0("> ! Hard-coded parameters:\n")
printAndLog(txt, logFile)
txt <- paste0("... caller = ",  caller, "\n")
printAndLog(txt, logFile)
txt <- paste0("... corMethod = ",  corMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("... distLimit = ",  distLimit, "\n")
printAndLog(txt, logFile)
txt <- paste0("... fitMeth = ",  fitMeth, "\n")
printAndLog(txt, logFile)

sameTADcol <- "darkorange1"
diffTADcol <- "darkslateblue"
mycols <- c("same TAD" ="darkorange1" , "diff. TAD"="darkslateblue")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

toprankingScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_topRanking")
source(paste0(toprankingScriptDir, "/", "get_topTADs.R"))

utilsDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")
source(file.path(utilsDir, "coreg_utils_ggscatterhist.R"))

################################################ DATA PREPARATION

if(buildTable) {
  
  cat(paste0("... load DIST data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_DIST_SORTNODUP/all_dist_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_dist_pairs)
  nrow(all_dist_pairs)
  all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
  all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
  ### UPDATE 30.06.2018
  stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
  
  cat(paste0("... load TAD data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_TAD_SORTNODUP/all_TAD_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_TAD_pairs)
  nrow(all_TAD_pairs)
  all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
  all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
  ### UPDATE 30.06.2018
  stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
  
  cat(paste0("... load COEXPR data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_COEXPR_SORTNODUP"),  paste0(curr_dataset, "_", corMethod), "coexprDT.Rdata"))
  cat(paste0(Sys.time(), "\n"))
  head(coexprDT)
  nrow(coexprDT)
  coexprDT$gene1 <- as.character(coexprDT$gene1)
  coexprDT$gene2 <- as.character(coexprDT$gene2)
  ### UPDATE 30.06.2018
  stopifnot(coexprDT$gene1 < coexprDT$gene2)
  
  #============================== RETRIEVE PIPELINE DATA FOR THIS DATASET
  dataset_pipDir <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)
  
  script0_name <- "0_prepGeneData"
  pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
  pipeline_geneList <- as.character(pipeline_geneList)
  
  
  dataset_dist_pair <- all_dist_pairs[all_dist_pairs$gene1 %in% pipeline_geneList & 
                                        all_dist_pairs$gene2 %in% pipeline_geneList,]
  
  dataset_dist_pairs_limit <- dataset_dist_pair[dataset_dist_pair$dist <= distLimit,]
  head(dataset_dist_pairs_limit)
  nrow(dataset_dist_pairs_limit)
  
  dataset_TAD_pairs <- all_TAD_pairs[all_TAD_pairs$gene1 %in% pipeline_geneList & 
                                       all_TAD_pairs$gene2 %in% pipeline_geneList,]

  # START MERGING DATA 
  
  cat(paste0("... merge DIST - TAD data\t", Sys.time(), "\t"))
  dataset_dist_TAD_DT <- left_join(dataset_dist_pairs_limit, dataset_TAD_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  dataset_dist_TAD_DT$sameTAD <- ifelse(is.na(dataset_dist_TAD_DT$region), 0, 1)

  cat(paste0("... merge COEXPR data\t", Sys.time(), "\t"))
  
  dataset_dist_TAD_coexpr_DT <- left_join(dataset_dist_TAD_DT, coexprDT, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  allData_dt <- dataset_dist_TAD_coexpr_DT
  allData_dt$region <- NULL
  allData_dt <- na.omit(allData_dt)
  
  outFile <-file.path(outFold, paste0("allData_dt.Rdata"))
  save(allData_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))    
}  else{
  outFile <-file.path(outFold, paste0("allData_dt.Rdata"))
  allData_dt <- eval(parse(text = load(outFile)))
  # load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3/TCGAcrc_msi_mss_hgnc/hgnc_family_allData_dt.Rdata")
}

stopifnot(nrow(allData_dt) > 0)

allData_dt$dist_kb <- allData_dt$dist/1000

sameTAD_DT <- allData_dt[allData_dt$sameTAD == 1,c("gene1", "gene2", "coexpr", "dist_kb")]
sameTAD_DT <- na.omit(sameTAD_DT)
sameTAD_DT <- sameTAD_DT[order(sameTAD_DT$dist_kb),]
sameTAD_DT$label <- "same TAD"

diffTAD_DT <- allData_dt[allData_dt$sameTAD == 0,c("gene1", "gene2",  "coexpr","dist_kb")]
diffTAD_DT <- na.omit(diffTAD_DT)
diffTAD_DT <- diffTAD_DT[order(diffTAD_DT$dist_kb),]
diffTAD_DT$label <- "diff. TAD"

sum_4_curves_DT <- rbind(rbind(sameTAD_DT, diffTAD_DT))

stopifnot(sum_4_curves_DT$dist_kb <= distLimit/1000)

sum_4_curves_DT$label <- factor(sum_4_curves_DT$label,
                                levels=names(mycols))

stopifnot(!any(is.na(sum_4_curves_DT$label)))


if(fitMeth == "loess") {
  
  sameTAD_DT <- allData_dt[allData_dt$sameTAD == 1, ]
  stopifnot(is.numeric(sameTAD_DT$coexpr[1]))
  stopifnot(is.numeric(sameTAD_DT$dist[1]))

  diffTAD_DT <- allData_dt[allData_dt$sameTAD == 0, ]
  stopifnot(is.numeric(diffTAD_DT$coexpr[1]))
  stopifnot(is.numeric(diffTAD_DT$dist[1]))
  
  my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
  my_xlab <- paste0("Distance between the 2 genes (kb)")
  my_sub <- paste0(curr_dataset)
  
  p_coexpr2curves <- my_ggscatterhist(
    sum_4_curves_DT,
    x = "dist_kb", 
    y = "coexpr",
    xlab=my_xlab,
    ylab=my_ylab,
    title = paste0("Coexpression along distance betw. genes"),
    subtitle=my_sub,
    legend.title="",
    point = FALSE,
    # rug=TRUE,
    add = fitMeth,
    color = "label", size = 3, alpha = 0.6,
    palette = mycols,
    x_margin.params = list(fill = "label", color = "black", size = 0.2),
    x_margin.ggtheme = theme_minimal(),
    x_margin.plot = "density",
    y_margin.params = list(fill = "label", color = "black", size = 0.2),
    y_margin.ggtheme = theme_minimal(),
    y_margin.plot = "boxplot",
    plot_xmargin = TRUE,
    plot_ymargin=TRUE,
    ymargin_as_xmargin = FALSE,
    # global_margins_cm = c(0, 0, 0.25, 0.25)
    global_margins_cm = c(0.25, 0.25, 0.25, 0.25)
  )
  outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_with_boxplot.", plotType))
  ggsave(p_coexpr2curves, height = myHeight*1.2, width = myWidth*1.2, filename=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  # PREDICT WITH ORIGINAL DISTANCE VALUES
  my_xlab <- paste0("Distance between the 2 genes (bp)")
  
  smooth_vals_sameTAD <- predict(loess(coexpr ~ dist, data = sameTAD_DT), sort(sameTAD_DT$dist))
  smooth_vals_diffTAD <- predict(loess(coexpr ~ dist, data = diffTAD_DT), sort(diffTAD_DT$dist))
  
  auc_diffTAD_obsDist <- auc(x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD)
  auc_sameTAD_obsDist <- auc(x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD)
  
  outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_originalDist.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(NULL,
       xlim = range(allData_dt$dist), 
       ylim = range(c(smooth_vals_sameTAD, smooth_vals_diffTAD)),
       # xlab="", 
       # ylab="",
       xlab=my_xlab, 
       ylab=my_ylab,
       main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
  mtext(text = "observed distance values", side = 3)
  lines( x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD, col = sameTADcol)
  lines( x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD, col = diffTADcol)
  legend("topright", 
         legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_obsDist, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_obsDist, 2))),  
         col = c(sameTADcol, diffTADcol),
         lty=1,
         bty = "n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  # PREDICT WITH DISTANCE VECTOR
  distVect <- seq(from=0, to = distLimit, length.out = nbrLoessPoints)
  smooth_vals_sameTAD_distVect <- predict(loess(coexpr ~ dist, data = sameTAD_DT), distVect)
  smooth_vals_diffTAD_distVect <- predict(loess(coexpr ~ dist, data = diffTAD_DT), distVect)

  auc_diffTAD_distVect <- auc(x = distVect, y = smooth_vals_diffTAD_distVect)
  auc_sameTAD_distVect <- auc(x = distVect, y = smooth_vals_sameTAD_distVect)

  
  diffTAD_mod <- loess(coexpr ~ dist, data = diffTAD_DT)
  outFile <- file.path(outFold, "diffTAD_mod.Rdata")
  save(diffTAD_mod, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  diffTAD_obsDist <- diffTAD_DT$dist
  outFile <- file.path(outFold, "diffTAD_obsDist.Rdata")
  save(diffTAD_obsDist, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFold, "auc_diffTAD_distVect.Rdata")
  save(auc_diffTAD_distVect, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))


  outFile <- file.path(outFold, "auc_sameTAD_distVect.Rdata")
  save(auc_sameTAD_distVect, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))


  outFile <- file.path(outFold, "smooth_vals_sameTAD_distVect.Rdata")
  save(smooth_vals_sameTAD_distVect, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))


  outFile <- file.path(outFold, "smooth_vals_diffTAD_distVect.Rdata")
  save(smooth_vals_diffTAD_distVect, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))




  outFile <- file.path(outFold, "distVect.Rdata")
  save(distVect, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))

  
  outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_vectDist.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(NULL,
       xlim = range(distVect), 
       ylim = range(c(na.omit(smooth_vals_sameTAD_distVect), na.omit(smooth_vals_diffTAD_distVect))),
       # xlab="", 
       # ylab="",
       xlab=my_xlab,
       ylab=my_ylab,
       main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
  mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
  lines( x = distVect, y = smooth_vals_sameTAD_distVect, col = sameTADcol)
  lines( x = distVect, y = smooth_vals_diffTAD_distVect, col = diffTADcol)
  legend("topright", 
         legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_distVect, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_distVect, 2))), 
         col = c(sameTADcol, diffTADcol),
         lty=1,
         bty = "n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  auc_values <- list(
    auc_diffTAD_distVect = auc_diffTAD_distVect,
    auc_sameTAD_distVect = auc_sameTAD_distVect,
    auc_ratio_same_over_diff_distVect = auc_sameTAD_distVect/auc_diffTAD_distVect,
    auc_diffTAD_obsDist = auc_diffTAD_obsDist,
    auc_sameTAD_obsDist = auc_sameTAD_obsDist,
    auc_ratio_same_over_diff_obsDist = auc_sameTAD_distVect/auc_diffTAD_obsDist
  )
  
  outFile <- file.path(outFold, paste0("auc_values.Rdata"))
  save(auc_values, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else{
  stop("only loess implemented yet\n")
  
}

  
######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
