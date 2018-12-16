startTime <- Sys.time()
cat(paste0("> Rscript scores_vs_other_variables_withBoxplots.R\n"))

#  Rscript scores_vs_other_variables_withBoxplots.R

library(foreach)
library(doMC)

source("analysis_utils.R")

signifThresh <- 0.05
signifVar <- "adj.P.Val"

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- TRUE

registerDoMC(ifelse(SSHFS, 2, 40))

printAndLog <- function(txt, logfile) {
  cat(txt)
  cat(txt, file = logfile, append=T)
}

all_datasets <- c("GSE87194_control_schi")

caller <- "TopDom"
script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11_name <- "11_runEmpPvalCombined"
script170_name <- "170_score_auc_pval_withShuffle"

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 10)
# myWidth <- ifelse(plotType == "png", 600, 10)
myWidth <- myHeight

outFold <- file.path("SCORES_VS_OTHER_VARIABLES_WITH_BOXPLOTS")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "score_vs_other_variables_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotAxisCex <- 1.2


aucCoexprDistFolder <- file.path(setDir, 
                                 "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUC_COEXPRDIST_SORTNODUP")
stopifnot(file.exists(aucCoexprDistFolder))

pipOutFolder <-  file.path(setDir, 
                           "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER")
stopifnot(file.exists(pipOutFolder))

settingFilesFolder <- file.path(setDir,
                                "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput"
                                )
stopifnot(file.exists(settingFilesFolder))

geneVarFile <- file.path(setDir,
                           "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18",
                            "GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata"
                           )

cat(paste0("!!! variance retrieved from: ", geneVarFile, "\n"))

stopifnot(file.exists(geneVarFile))

all_datasets <- list.files(pipOutFolder)
# head(all_datasets)
stopifnot(length(all_datasets) >  0 )
cat(paste0("... found ", length(all_datasets), " datasets\n"))

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

cancerDS <- score_DT$dataset[score_DT$process_short == "cancer"]
noTCGA_cancerDS <- cancerDS[!grepl("^TCGA", cancerDS)]
TCGA_cancerDS <- cancerDS[grepl("TCGA", cancerDS)]
no_cancerDS <- score_DT$dataset[score_DT$process_short != "cancer"]


# all_datasets = all_datasets[1]
# all_datasets = all_datasets[1:5]

curr_dataset = "TCGAcrc_msi_mss"

# return NULL if file not found ??? [if yes -> build table skipping missing files, otherwise raise error and stop]
returnNull <-  FALSE

txt <- paste0("... geneVarFile\t=\t", geneVarFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... returnNull\t=\t", as.character(returnNull), "\n")
printAndLog(txt, logFile)
txt <- paste0("... found # noTCGA_cancerDS\t=\t", length(noTCGA_cancerDS) , "\n" )
printAndLog(txt, logFile)
txt <- paste0("... found # TCGA_cancerDS\t=\t", length(TCGA_cancerDS) , "\n" )
printAndLog(txt, logFile)
txt <- paste0("... found # no_cancerDS\t=\t", length(no_cancerDS) , "\n" )
printAndLog(txt, logFile)

if(buildTable) {
  
  geneVarDT <- eval(parse(text = load(geneVarFile)))
  
  datasets_variables_DT <- foreach(curr_dataset = all_datasets, .combine='rbind') %dopar% {
    
    # txt <- paste0("> START ", curr_dataset, "\n")
    # printAndLog(txt, logFile)
    
    cat("... load auc coexpr dist\n")
    aucCoexprDistFile <- file.path(aucCoexprDistFolder,
                                   curr_dataset, "auc_values.Rdata")
    if(returnNull) { if(!file.exists(aucCoexprDistFile)) return(NULL) } else { stopifnot(file.exists(aucCoexprDistFile)) }
    all_coexprDistAUC <- eval(parse(text = load(aucCoexprDistFile)))
    coexprDistAUC <- all_coexprDistAUC[["auc_ratio_same_over_diff_distVect"]]
    cat("coexprDistAUC = ", coexprDistAUC, "\n")
    stopifnot(is.numeric(coexprDistAUC))
    
    cat("... load auc FCC \n")
    aucFCC_file <- file.path(pipOutFolder, curr_dataset, script170_name, "allratio_auc_pval.Rdata")

    if(returnNull) { if(!file.exists(aucFCC_file)) return(NULL) } else { stopifnot(file.exists(aucFCC_file)) }

    all_fccAUC <- eval(parse(text = load(aucFCC_file)))
    
    fccAUC <- all_fccAUC["prodSignedRatio_auc_permGenes"]
    stopifnot(is.numeric(fccAUC))
    
    cat("... load pipeline_geneList\n")
    geneFile <- file.path(pipOutFolder, 
                          curr_dataset, 
                          script0_name,
                          "pipeline_geneList.Rdata")
    # if(!file.exists(geneFile)) return(NULL)
    pipeline_geneList <- eval(parse(text = load(geneFile)))
    
    tmp_g2t <- gene2tad_DT[gene2tad_DT$entrezID %in% pipeline_geneList,]
    stopifnot(nrow(tmp_g2t) > 0)
    nGenesByTAD_dt <- data.frame(region = as.character(names(table(tmp_g2t$region))),
                                 nGenes = as.numeric(table(tmp_g2t$region)),
                                 stringsAsFactors = FALSE
                                 )
    meanGenesByTAD <- mean(nGenesByTAD_dt$nGenes)
    medianGenesByTAD <- median(nGenesByTAD_dt$nGenes)
    
    
    cat("... load limma DT\n")
    limmaFile <- file.path(pipOutFolder, 
                           curr_dataset, 
                           script1_name,
                           "DE_topTable.Rdata")

    if(returnNull) { if(!file.exists(limmaFile)) return(NULL) } else { stopifnot(file.exists(limmaFile)) }

    limmaDT <- eval(parse(text = load(limmaFile)))
    
    stopifnot(limmaDT$genes == rownames(limmaDT))
    stopifnot(names(pipeline_geneList) %in% limmaDT$genes)
    
    limmaDT <- limmaDT[limmaDT$genes %in% names(pipeline_geneList),]
    stopifnot(nrow(limmaDT) > 0)

    cat("... load TAD\n")
    tadpvalFile <-  file.path(pipOutFolder, curr_dataset, script11_name, "emp_pval_combined.Rdata")

    if(returnNull) { if(!file.exists(tadpvalFile)) return(NULL) } else { stopifnot(file.exists(tadpvalFile)) }

    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    
    
    cat("... load setting files\n")
    settingF <- file.path(settingFilesFolder,
                          paste0("run_settings_", curr_dataset, ".R"))

    if(returnNull) { if(!file.exists(settingF)) return(NULL) } else { stopifnot(file.exists(settingF)) }

    source(settingF)
    sample1_file <- file.path(setDir, sample1_file)
    if(returnNull) { if(!file.exists(sample1_file)) return(NULL) } else { stopifnot(file.exists(sample1_file)) }
    s1 <- eval(parse(text = load(sample1_file)))
    
    sample2_file <- file.path(setDir, sample2_file)
    if(returnNull) { if(!file.exists(sample2_file)) return(NULL) } else { stopifnot(file.exists(sample2_file)) }
    s2 <- eval(parse(text = load(sample2_file)))
    
    nGenes <- length(pipeline_geneList)
    nSignifGenes <- sum(limmaDT[,signifVar] <= signifThresh)
    
    nTADs <- length(tad_pval)
    nSignifTADs <- sum(tad_pval <= signifThresh)
    
    nSamp1 <- length(s1)
    nSamp2 <- length(s2)
    
    stopifnot(curr_dataset %in% geneVarDT$dataset)
    meanMostVar <- geneVarDT$meanMostVar[geneVarDT$dataset == curr_dataset]

    data.frame(
      dataset=curr_dataset,
      coexprDistAUC = coexprDistAUC,
      fccAUC = fccAUC,
      nSamp1 = nSamp1,
      nSamp2 = nSamp2,
      nAllSamp = nSamp1+nSamp2,
      ratioSamp = nSamp1/nSamp2,
      nGenes = nGenes,
      nSignifGenes = nSignifGenes,
      nTADs = nTADs,
      nSignifTADs=nSignifTADs,
      meanGenesByTAD = meanGenesByTAD,
      medianGenesByTAD = medianGenesByTAD,
      meanMostVar = meanMostVar,
      stringsAsFactors = FALSE
      
    )
    
  }
  rownames(datasets_variables_DT) <- NULL
  stopifnot(nrow(datasets_variables_DT) == length(all_datasets))
  outFile <- file.path(outFold, "datasets_variables_DT.Rdata")
  save(datasets_variables_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFold, "datasets_variables_DT.Rdata")
  #outFile <- "SCORES_VS_OTHER_VARIABLES_WITH_BOXPLOTS/datasets_variables_DT.Rdata"
  stopifnot(file.exists(outFile))
  datasets_variables_DT <- eval(parse(text = load(outFile)))
}

# put it before adding color and type columns...
all_vars <- colnames(datasets_variables_DT)[! colnames(datasets_variables_DT) %in% c("dataset", "nSamp1", "nSamp2")]

datasets_variables_DT$color <- ifelse(datasets_variables_DT$dataset %in% noTCGA_cancerDS, dataset_colors["noTCGA_cancerDS"],
                                      ifelse(datasets_variables_DT$dataset %in% TCGA_cancerDS, dataset_colors["TCGA_cancerDS"],
                                             ifelse(datasets_variables_DT$dataset %in% no_cancerDS, dataset_colors["no_cancerDS"], NA)))
stopifnot(!is.na(datasets_variables_DT$color))


datasets_variables_DT$dsType <- ifelse(datasets_variables_DT$dataset %in% noTCGA_cancerDS, "cancer_notTCGA",
                                   ifelse(datasets_variables_DT$dataset %in% TCGA_cancerDS, "TCGA",
                                          ifelse(datasets_variables_DT$dataset %in% no_cancerDS, "not_cancer", NA)))

stopifnot(!is.na(datasets_variables_DT$dsType))

datasets_variables_DT$dsType <- factor(datasets_variables_DT$dsType, levels = c("not_cancer", "cancer_notTCGA","TCGA"))
stopifnot(!is.na(datasets_variables_DT$dsType))
datasets_variables_DT <- datasets_variables_DT[order(as.numeric(datasets_variables_DT$dsType)),]

# stop("--ok\n")

######################################################################################
###################################################################################### # draw the boxplots
######################################################################################
# dataset coexprDistAUC   fccAUC nSamp1 nSamp2 nAllSamp ratioSamp nGenes nSignifGenes nTADs nSignifTADs meanGenesByTAD medianGenesByTAD meanMostVar
# GSE101521_control_mdd      1.130653 1.136174     29     21       50  1.380952  14869            0  2459         402       6.046767                5   0.5967525
# GSE102073_stic_nostic      1.187492 1.506043     43     42       85  1.023810  14747            0  2436        1946       6.053777                5   5.0150059

curr_var  = "nAllSamp"

datasets_variables_DT$dsTypeLabel <- gsub("_", "\n", as.character(datasets_variables_DT$dsType))
datasets_variables_DT$dsTypeLabel <- factor(datasets_variables_DT$dsTypeLabel, levels = gsub("_", "\n", levels(datasets_variables_DT$dsType)))
stopifnot(!is.na(datasets_variables_DT$dsTypeLabel))


for(curr_var in all_vars) {
  
  outFile <- file.path(outFold, paste0(curr_var, "_", "cancer_TCGA_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(as.formula(paste0(curr_var, " ~  dsTypeLabel")), 
          data = datasets_variables_DT,
          boxfill = unique(datasets_variables_DT$color),
          main = paste0(curr_var),
          ylab = paste0(curr_var),
          las = 2,
          cex.lab = plotAxisCex, cex.axis = plotAxisCex
  )  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


######################################################################################
######################################################################################

stopifnot(!duplicated(datasets_variables_DT$dataset))
rownames(datasets_variables_DT) <- datasets_variables_DT$dataset
datasets_variables_DT$dataset <- NULL

other_vars <- all_vars

titNames <- c(
  coexprDistAUC="AUC ratio - pairwise coexpr.",
  fccAUC = "AUC ratio - FCC",
  nSamp1 = "# samples (cond1)",
  nSamp2 = "# samples (cond2)",
  nAllSamp = "# samples (all)",
  ratioSamp = "ratio cond1/cond2 samples",
  nGenes = "# genes",
  nSignifGenes = paste0("# signif. DE genes (", signifVar, " <= ", signifThresh, ")"),
  nTADs = "# TADs",
  nSignifTADs = paste0("# signif. DA TADs (", signifVar, " <= ", signifThresh, ")"),
  meanGenesByTAD = "mean # genes/TAD",
  medianGenesByTAD = "median # genes/TAD",
  meanMostVar = "mean var. of most var. genes"
)

offSets <- c(
  coexprDistAUC=0.03,
  fccAUC = 0.03,
  nSamp1 = 150,
  nSamp2 = 150,
  nAllSamp = 150,
  ratioSamp = 5,
  nGenes = 1000,
  nSignifGenes = 1000,
  nTADs = 100,
  nSignifTADs = 500,
  meanGenesByTAD = 2,
  medianGenesByTAD = 2,
  meanMostVar = 0.1
)


# stopifnot(colnames(datasets_variables_DT) %in% names(titNames))

mySub <- paste0("all datasets (n=", nrow(datasets_variables_DT), ")")

stopifnot(rownames(datasets_variables_DT) %in% names(dataset_proc_colors))
curr_colors <- dataset_proc_colors[rownames(datasets_variables_DT)]

for(ref_var in c("coexprDistAUC", "fccAUC")) {
  for(curr_var in other_vars ) {
    
    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, ".", plotType))
    
    myx <- datasets_variables_DT[,curr_var]
    myy <- datasets_variables_DT[,ref_var]
    
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(varX = myx, varY=myy, 
                     mylabels = as.character(rownames(datasets_variables_DT)), 
                     withLab=F,
                     main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                     xlab = paste0(titNames[curr_var]),
                     ylab = paste0(titNames[ref_var]),
                     xlim = range(myx) + c(-offSets[curr_var], offSets[curr_var]),
                     ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                     cex.lab = plotAxisCex, cex.axis = plotAxisCex
                     )
    text(x = myx, y = myy,
         labels =  as.character(rownames(datasets_variables_DT)), 
         pch=16,
         col = curr_colors,
         bty="n",
         cex=0.7)
    mtext(side=3, text = mySub)
    legend("bottomright",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    if(grepl("^n.*Genes", curr_var) | grepl("TADs", curr_var)) {
      myx <- datasets_variables_DT[,curr_var]
      myy <- datasets_variables_DT[,ref_var]
      
      myx <- log10(myx+1)
      
      if(grepl("^n.*Genes", curr_var)) {
      xoffset <- ifelse(curr_var == "nGenes", 0.05,
                        ifelse(curr_var == "nSignifGenes", 0.35, stop("")))
      } else if(grepl("TADs", curr_var)) {
        xoffset <- ifelse(curr_var == "nTADs", 0.05,
                          ifelse(curr_var == "nSignifTADs", 0.35, stop("")))
      } else {stop("")}
      outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, "_log10.", plotType))
      
      do.call(plotType, list(outFile, height = myHeight, width = myWidth))
      my_plot_function(varX = myx, varY=myy, 
                       mylabels = as.character(rownames(datasets_variables_DT)), 
                       withLab=F,
                       main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                       xlab = paste0(titNames[curr_var], "(log10[#+1])" ),
                       ylab = paste0(titNames[ref_var]),
                       xlim = range(myx) + c(-xoffset, xoffset),
                       ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                       cex.lab = plotAxisCex, cex.axis = plotAxisCex
      )
      text(x = myx, y = myy,
           labels =  as.character(rownames(datasets_variables_DT)), 
           pch=16,
           col = curr_colors,
           bty="n",
           cex=0.7)
      mtext(side=3, text = mySub)
      legend("bottomright",
             legend=names(my_colors),
             lty=1,
             col = my_colors,
             lwd = 5,
             bty="n",
             cex = 0.7)
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
    }
    
    
    
    
  } # end iterating over variables
} # end iterating FCC - coexprdist


# just to check: # DE genes and meanMostVar

#all_curr_vars <- c("nSignifGenes", "nSignifTADs")
#ref_var <- "meanMostVar" 

all_cmps <- list(
c(curr_var = "nSignifGenes", ref_var= "meanMostVar"),
c(curr_var = "nSignifTADs", ref_var= "meanMostVar"),
c(curr_var = "nSignifGenes", ref_var= "nSignifTADs")
)

for(i in seq_along(all_cmps)){

    curr_var = all_cmps[[i]]["curr_var"]
    ref_var = all_cmps[[i]]["ref_var"]

    myx <- log10(datasets_variables_DT[, curr_var]+1)
    myy <- datasets_variables_DT[, ref_var]

    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, "_log10.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(varX = myx, varY=myy, 
                     mylabels = as.character(rownames(datasets_variables_DT)), 
                     withLab=F,
                     main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                     xlab = paste0(titNames[curr_var], "(log10[#+1])" ),
                     ylab = paste0(titNames[ref_var]),
                     xlim = range(myx) + c(-1,1),# + c(-log10(xoffset[curr_var]), log10(xoffset[curr_var])),
                     ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                     cex.lab = plotAxisCex, cex.axis = plotAxisCex
    )
    text(x = myx, y = myy,
         labels =  as.character(rownames(datasets_variables_DT)), 
         pch=16,
         col = curr_colors,
         bty="n",
         cex=0.7)
    mtext(side=3, text = mySub)
    legend("bottomright",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))

    myx <- datasets_variables_DT[, curr_var]
    myy <- datasets_variables_DT[, ref_var]

    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(varX = myx, varY=myy, 
                     mylabels = as.character(rownames(datasets_variables_DT)), 
                     withLab=F,
                     main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                     xlab = paste0(titNames[curr_var], "" ),
                     ylab = paste0(titNames[ref_var]),
                     xlim = range(myx) + c(-1000,1000),#+ c(-xoffset[curr_var], xoffset[curr_var]),
                     ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                     cex.lab = plotAxisCex, cex.axis = plotAxisCex
    )
    text(x = myx, y = myy,
         labels =  as.character(rownames(datasets_variables_DT)), 
         pch=16,
         col = curr_colors,
         bty="n",
         cex=0.7)
    mtext(side=3, text = mySub)
    legend("bottomright",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
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
