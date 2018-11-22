#!/usr/bin/Rscript

cat(paste0("> START ", "wave_plot_aucFCC.R",  "\n"))


startTime <- Sys.time()

# Rscript wave_plot_aucFCC.R 
# Rscript wave_plot_aucFCC.R ../TAD_DE_pipeline/SETTING_FILES_cleanDE/run_settings_TCGAcrc_msi_mss.R
# =>> NEED TO ADD curr_ds IN OUTPUT NAME !!!

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

nCpu <- ifelse(SSHFS, 2, 30)

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

### SET OUTPUT
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)

# if plotSeparated == TRUE => 1 plot per ratio, otherwise all on the same figure (#x2 plots)
plotSeparated <- F

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  stopifnot(length(args) == 1)
  all_setting_files <- args
} else {
  settingFilesFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")
  stopifnot(file.exists(settingFilesFold))
  all_setting_files <- list.files(settingFilesFold, full.names = TRUE, pattern = paste0("run_settings_.+\\.R$"))
} 
stopifnot(file.exists(all_setting_files))

cat(paste0("... will analyse # datasets = ", length(all_setting_files), "\n"))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
pipRunDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8c_runAllDown"

source("analysis_utils.R")
source(file.path(pipRunDir, "main_settings.R"))
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

# create the directories
outFold <- paste0("WAVE_PLOT_AUCFCC")
system(paste0("mkdir -p ", outFold))


logFile <- file.path(outFold, "wave_plot_aucfcc_logfile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""


# 14f:   leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))
# 14f2:  leg_xpos <- observed_auc

allDown_limited <- "prodSignedRatio"

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

init_gene2tadDT <- gene2tadDT

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# all_setting_files <- all_setting_files[1:2]

all_plotVect <- list()


for(settingF in all_setting_files) {
  
  
  
  source(settingF)
  
  curr_ds <- gsub("run_settings_(.+)\\.R", "\\1", basename(settingF))
  
  pipOutFold <- file.path(pipRunDir, pipOutFold)
  
  
  # UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
  initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
  geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
  txt <- paste0("> Start with # genes: ", length(geneList), "/", length(initList), "\n")
  printAndLog(txt, logFile)
  
  stopifnot(!any(duplicated(names(geneList))))
  
  gene2tadDT <- init_gene2tadDT[init_gene2tadDT$entrezID %in% geneList,]
  geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))
  

  ###############
  ##### retrieve the direction of up/down
  ###############
  # retrieve the direction of up/down
  DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
  DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
  exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
  # samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
  # samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
  samp1 <- eval(parse(text=load(paste0(sample1_file))))
  samp2 <- eval(parse(text=load(paste0(sample2_file))))
  DE_topTable <- DE_topTable[DE_topTable$genes %in% names(DE_geneList),]
  stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
  stopifnot(!any(duplicated(names(DE_geneList))))
  stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
  stopifnot(all(samp1 %in% colnames(exprDT)))
  stopifnot(all(samp2 %in% colnames(exprDT)))
  maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
  stopifnot(maxDownGene %in% rownames(exprDT))
  mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
  mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)
  
  if(mean_expr1 > mean_expr2) {
    subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
  } else{
    subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
  }
  
#  if(! plotSeparated) {
#    nColPlot <- 2
#    nRowPlot <- length(allDown_limited)*2/nColPlot
#    outFile <- file.path(outFold, paste0(curr_ds, "_allRatios_cumsum_obs_permut.", plotType))
#    do.call(plotType, list(outFile, height=myHeight*nRowPlot, width=myWidth*nColPlot))
#    par(mfrow=c(nRowPlot, nColPlot))
#  }


  outFile <- file.path(outFold, paste0(curr_ds, "_allRatios_cumsum_obs_permut.", plotType))
  do.call(plotType, list(outFile, height=myHeight*1, width=myWidth*1))

  
  ################################****************************************************************************************
  ####################################################### ITERATE OVER RATIOS TO PLOT
  ################################****************************************************************************************
  
  for(curr_ratio_type in allDown_limited) {
    cat(paste0("*** START ", curr_ratio_type, "\n"))
    
    obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
    permut_currDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
    
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
    
    interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )
    
    ############################################################
    ############################################################ # filter the TADs and sort
    ############################################################
    filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)
    
    filter_permut_currDown_unsort <- permut_currDown[interReg,]
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
    
    filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
    rownames(filter_permut_currDown) <- NULL
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
    
    # FOR ratioDown => plot ratioConcord, departure from 0.5
    if(curr_ratio_type == "ratioDown") {
      my_stat_curr_ratio <- "ratioDown_Concord"
      departureFromValue <- 0.5
      # => Concord, departure 0.5
      # Transform ratioDown -> ratioConcord
      # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
      # transf. ratioDown -> ratioConcord
      filter_obs_curr_down_half <- abs(filter_obs_curr_down - 0.5) + 0.5
      filter_permut_currDown_half <- abs(filter_permut_currDown - 0.5) + 0.5

      my_main <- paste0(curr_ratio_type, ": cumsum departure from ", departureFromValue)
      my_ylab <- paste0("cumsum(abs(", curr_ratio_type, " - ", departureFromValue,"))")
      my_xlab <- paste0("regions ranked by decreasing ", curr_ratio_type)


    } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
      my_stat_curr_ratio <- paste0(curr_ratio_type, "_Concord")
      departureFromValue <- 0.5
      # => Concord, departure 0.5
      # Transform rescWeightedQQ -> rescWeightedQQConcord
      # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
      # transf. ratioDown -> ratioConcord
      filter_obs_curr_down_half <- abs(filter_obs_curr_down - 0.5) + 0.5
      filter_permut_currDown_half <- abs(filter_permut_currDown - 0.5) + 0.5

      my_main <- paste0(curr_ratio_type, ": cumsum departure from ", departureFromValue)
      my_ylab <- paste0("cumsum(abs(", curr_ratio_type, " - ", departureFromValue,"))")
      my_xlab <- paste0("regions ranked by decreasing ", curr_ratio_type)


    } else if(curr_ratio_type == "prodSignedRatio") {
      my_stat_curr_ratio <- "prodSignedRatio"
      departureFromValue <- 0
      # => raw (departure 0)
      # prodSignedRatio -> does not need to be transformed
      filter_obs_curr_down_half <- filter_obs_curr_down
      filter_permut_currDown_half <- filter_permut_currDown
      curr_ratio_typeTit <- "FCC score"
      my_main <- paste0(curr_ds, ": ", curr_ratio_typeTit, " cumsum")
      my_ylab <- paste0("cumsum(abs(", curr_ratio_typeTit, "))")
      my_xlab <- paste0("TADs ranked by decreasing ", curr_ratio_typeTit)

    }  
    # PLOT THE 1ST PLOT
#    if(plotSeparated) {
#      outFile <- file.path(outFold, paste0(curr_ds, "_", curr_ratio_type, "_departure05_cumsum_obs_permut.", plotType))
#      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#    }
    plotVect <- plot_cumsumDiff05_returnVect(filter_obs_curr_down_half, 
                      filter_permut_currDown_half, 
                      my_stat = my_stat_curr_ratio,
                      departureValue = departureFromValue, drawline=TRUE,
                      main = my_main, xlab = my_xlab, ylab = my_ylab)
    mtext(subtitDir, font=3)

#    if(plotSeparated) {
#      foo <- dev.off()
#      cat(paste0("... written: ", outFile, "\n"))
#    }
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    all_plotVect[[paste0(curr_ds)]] <- plotVect
    
    # PLOT THE 2ND PLOT
#    if(plotSeparated){
#      outFile <- file.path(outFold, paste0(curr_ds, "_", curr_ratio_type, "_departure05_cumsum_obs_permut_AUC.", plotType))
#      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#    }
#    plot_cumsumDiff05_AUC2(filter_obs_curr_down_half, 
#                           filter_permut_currDown_half, 
#                           my_stat = my_stat_curr_ratio,
#                           departureValue = departureFromValue)
#    mtext(subtitDir, font=3)
#    if(plotSeparated){
#      foo <- dev.off()
#      cat(paste0("... written: ", outFile, "\n"))
#    }
  }
#  
#  if(!plotSeparated){
#    foo <- dev.off()
#    cat(paste0("... written: ", outFile, "\n"))
#  }
  
  
  
  
}

outFile <- file.path(outFold, "all_plotVect.Rdata")
save(all_plotVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
######################################################################################
######################################################################################

cat(paste0("*** DONE"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)




