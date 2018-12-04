# Rscript gene_diffVariance.R <exprType>
# Rscript gene_diffVariance.R fpkm
# Rscript gene_diffVariance.R log2fpkm
# Rscript gene_diffVariance.R NBvar
# Rscript gene_diffVariance.R voom
cat("> START: gene_diffVariance.R\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

library(foreach)
library(doMC)
library(tools)
library(DESeq2)

startTime <- Sys.time()

pointPch <- 16
pointCex <- 1
cexAxis <- 1
cexLab <- 1

rangeOffset <- 0.15

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)


registerDoMC(ifelse(SSHFS,2,40))

args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) > 0)
exprType <- "log2fpkm"
nTopLast <- 1000
if(!is.na(args[1]))
  exprType <- args[1]
if(!is.na(args[2]))
  nTopLast <- as.numeric(args[2])
stopifnot(is.numeric(nTopLast))
stopifnot(exprType %in% c("fpkm", "log2fpkm", "NBvar", "voom"))

exprTypeName <- ifelse(exprType == "fpkm", "FPKM", 
				  ifelse(exprType == "log2fpkm", "log2 FPKM",
					ifelse(exprType == "NBvar", "negative binomial var.",
						ifelse(exprType == "voom", "voom transf.", NA))))
stopifnot(!is.na(exprTypeName))

cat("... START with exprType =\t", exprType, "\n")
cat("... START with nTopLast =\t", nTopLast, "\n")

buildTable <- F

# registerDoMC(ifelse(SSHFS, 2, 20))

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))

settingFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

pipMainFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

dsFold <- file.path(pipMainFolder, "OUTPUT_FOLDER")

all_setting_files <- list.files(settingFolder, full.names=T)
all_setting_files <- all_setting_files[grep(".R$", all_setting_files)]
stopifnot(length(all_setting_files) > 0)

outFold <- file.path(paste0("GENE_DIFFVARIANCE"), toupper(exprType))
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

# all_setting_files <- all_setting_files[1:5]

if(buildTable) {
  all_ds_geneVarDT <- foreach(ds_file = all_setting_files, .combine="rbind") %dopar% {
  
    stopifnot(file.exists(ds_file))
    cat("... source settingFile", basename(ds_file), "\n")
    source(ds_file)
  
    curr_ds <- basename(pipOutFold)
    cat("... START", curr_ds, "\n")
  
    ds_pipFolder <- file.path(pipMainFolder, pipOutFold)
    cat(ds_pipFolder,"\n")
    stopifnot(file.exists(ds_pipFolder))
    
    cat("... load samp1\n")
    samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
    cat("... load samp2\n")
    samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))
    
    cat("... load geneList\n")
    geneList <- eval(parse(text = load(
      file.path(ds_pipFolder, "0_prepGeneData", "pipeline_geneList.Rdata")
    )))
    
    cat("... load exprDT\n")
    
    if(exprType == "fpkm" | exprType == "log2fpkm" | exprType == "NBvar") {
      exprDT <- eval(parse(text = load(
        file.path(ds_pipFolder, "0_prepGeneData", paste0("rna_", "fpkm", "DT.Rdata"))
      )))
    } else if(exprType == "voom") {
      exprDT <- eval(parse(text = load(
        file.path(ds_pipFolder, "1_runGeneDE", paste0(exprType, "_lmFitInputDT.Rdata"))
      )))
    } else{
      stop(paste0("!!! unimplemented exprType", exprType, "\n"))
    }
    
    stopifnot(samp1 %in% colnames(exprDT))
    stopifnot(samp2 %in% colnames(exprDT))
    stopifnot(names(geneList) %in% rownames(exprDT))
    
    curr_exprDT <- exprDT[rownames(exprDT) %in% names(geneList), c(samp1,samp2)]  
    stopifnot(is.numeric(curr_exprDT[1,1]))
    
    if(exprType == "log2fpkm") {
      curr_exprDT <- log2(curr_exprDT + 1)
      stopifnot(is.numeric(curr_exprDT[1,1]))
    } else if(exprType == "NBvar") {
        cat("!!! WARNING to use \"DESeq2:::varianceStabilizingTransformation\" the expression data are rounded !!!\n")
        curr_exprDT <- try(DESeq2:::varianceStabilizingTransformation(round(as.matrix(curr_exprDT))))
        if(class(curr_exprDT) == "try-error"){
          return(data.frame(
              nTopLast = nTopLast,
              data_type = exprType,
              dataset = curr_ds,
              cond1 = cond1,
              cond2 = cond2,
              meanDiffVar = NA,
              meanDiffVarMostVarAll = NA,
              meanDiffVarMostVarCond1 =NA,
              meanDiffVarMostVarCond2 = NA,
              stringsAsFactors = FALSE
          ))
        }
          
        stopifnot(is.numeric(curr_exprDT[1,1]))
      }
    
    curr_exprDT_cond1 <- curr_exprDT[, samp1]
    
    stopifnot(ncol(curr_exprDT_cond1) == length(samp1))
    stopifnot(nrow(curr_exprDT_cond1) == nrow(curr_exprDT))
    
    curr_exprDT_cond2 <- curr_exprDT[, samp2]
    
    stopifnot(ncol(curr_exprDT_cond2) == length(samp2))
    stopifnot(nrow(curr_exprDT_cond2) == nrow(curr_exprDT))
    
    geneVarAll <- apply(curr_exprDT, 1,  var, na.rm=T)
    stopifnot(length(geneVarAll) >= nTopLast)
    
    geneVarAll <- sort(geneVarAll, decreasing = TRUE)
    stopifnot(length(geneVarAll) >= nTopLast)
    mostVariantAll <- names(geneVarAll[1:nTopLast])
    
    geneVarCond1 <- apply(curr_exprDT_cond1, 1,  var, na.rm=T)
    geneVarCond1_init <- geneVarCond1
    stopifnot(length(geneVarCond1) >= nTopLast)
    
    geneVarCond1 <- sort(geneVarCond1, decreasing = TRUE)
    stopifnot(length(geneVarCond1) >= nTopLast)
    mostVariantCond1 <- names(geneVarCond1[1:nTopLast])
    
    geneVarCond2 <- apply(curr_exprDT_cond2, 1,  var, na.rm=T)
    geneVarCond2_init <- geneVarCond2
    stopifnot(length(geneVarCond2) >= nTopLast)
    
    geneVarCond2 <- sort(geneVarCond2, decreasing = TRUE)
    stopifnot(length(geneVarCond2) >= nTopLast)
    mostVariantCond2 <- names(geneVarCond2[1:nTopLast])
    
    stopifnot(length(mostVariantAll) == nTopLast)
    stopifnot(length(mostVariantCond1) == nTopLast)
    stopifnot(length(mostVariantCond2) == nTopLast)
    
    stopifnot(names(geneVarCond1_init) == names(geneVarCond2_init))
    meanDiffVar <- mean(abs(geneVarCond2_init-geneVarCond1_init), na.rm=T)
    
    stopifnot(mostVariantAll %in% names(geneVarCond1_init))
    stopifnot(mostVariantAll %in% names(geneVarCond2_init))
    meanDiffVarMostVarAll <- mean(abs(geneVarCond2_init[mostVariantAll] - geneVarCond1_init[mostVariantAll]), na.rm=T)
      
    stopifnot(mostVariantCond1 %in% names(geneVarCond1_init))
    stopifnot(mostVariantCond1 %in% names(geneVarCond2_init))
    
    meanDiffVarMostVarCond1 <- mean(abs(geneVarCond2_init[mostVariantCond1] - geneVarCond1_init[mostVariantCond1]), na.rm=T)
      
    stopifnot(mostVariantCond2 %in% names(geneVarCond1_init))
    stopifnot(mostVariantCond2 %in% names(geneVarCond2_init))
    
    meanDiffVarMostVarCond2 <- mean(abs(geneVarCond2_init[mostVariantCond2] - geneVarCond1_init[mostVariantCond2]), na.rm=T)
      
    
    data.frame(
      nTopLast = nTopLast,
      data_type = exprType,
      dataset = curr_ds,
      cond1 = cond1,
      cond2 = cond2,
      
      meanDiffVar = meanDiffVar,
      meanDiffVarMostVarAll = meanDiffVarMostVarAll,
      meanDiffVarMostVarCond1 =meanDiffVarMostVarCond1,
      meanDiffVarMostVarCond2 = meanDiffVarMostVarCond2,
        
        
        
      stringsAsFactors = FALSE
    )
    
  }
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  save(all_ds_geneVarDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  all_ds_geneVarDT <- eval(parse(text = load(outFile)))
}

# stop("--ok--\n")



add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]  + R2shiftY
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         font=3,
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}



########################################################################################## RETRIEVE FCC AND COEXPRDIST
##########################################################################################

all_ds <- unique(all_ds_geneVarDT$dataset)

aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
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
names(aucFCC) <- all_ds

stopifnot(names(aucFCC) %in% names(dataset_proc_colors) )
curr_colors <- dataset_proc_colors[names(aucFCC)]

aucCoexprDist <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
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
names(aucCoexprDist) <- all_ds


##########################################################################################
##########################################################################################

if(exprType == "NBvar") {
  all_ds_geneVarDT <- all_ds_geneVarDT[!is.na(all_ds_geneVarDT$meanMostVar) & !is.na(all_ds_geneVarDT$meanMostVar),]
}


rangeOffset_auc <- 0.2
rangeOffset_var <- 1

####### scatter plot most var, least var, and ratio

vars_to_plot <- colnames(all_ds_geneVarDT)[! colnames(all_ds_geneVarDT) %in% c("nTopLast", "data_type", "dataset", "cond1", "cond2")]

all_auc <- c("FCC", "CoexprDist")


mySubTit <- paste0("all datasets (n=", nrow(all_ds_geneVarDT), ")")

auc_type = "FCC"

for(auc_type in all_auc) {
  
  myylab <-  paste0("% AUC increase - ", auc_type)

  plot_var = "meanDiffVar"
  for(plot_var in vars_to_plot) {
    
    for(transf in c("", "log10")) {
      
    
    cat("... start", auc_type, "\n")
    
    curr_auc <- eval(parse(text = paste0("auc", auc_type)))
    stopifnot(all_ds_geneVarDT$dataset %in% names(curr_auc))
    all_ds_geneVarDT[, auc_type] <- curr_auc[all_ds_geneVarDT$dataset]
    
    myx <- all_ds_geneVarDT[, plot_var]
    myy <- all_ds_geneVarDT[, auc_type]
    
    myxlab <- paste0(plot_var)
    
    if(grepl("MostVar", plot_var)) {
      myTit <- paste0(auc_type, ": ", plot_var, " (nTop=", nTopLast, "; exprType=", exprType, ")")  
    } else {
      myTit <- paste0(auc_type, ": ", plot_var, " (exprType=", exprType, ")")  
    }
    
    if(transf == "log10"){
      outFile <- file.path(outFold, paste0(auc_type, "auc_vs_", plot_var,  "_", transf, ".", plotType))  
      myx <- log10(myx)
      myxlab <- paste0(plot_var, " (log10)")
    } else{
      outFile <- file.path(outFold, paste0(auc_type, "auc_vs_", plot_var,  ".", plotType))
    }
    
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(
      x =  myx,
      y =  myy,
      xlim =  range(na.omit(myx) + c(-rangeOffset_var,rangeOffset_var)),
      ylim =  range(na.omit(myy) + c(-rangeOffset_auc,rangeOffset_auc)),
      pch = 16, cex = 0.7,
      col = curr_colors,
      ylab = myylab,
      xlab = myxlab,
      main = myTit
    )
    text(x = myx,
         y = myy,
         labels = all_ds_geneVarDT$dataset,
         col = curr_colors,
         pos=3, cex = 0.7)
    add_legend_narm(x = myx, 
                    y = myy,
                    mypos="bottomright", mymet="Pearson")  
    
    mtext(text = mySubTit, side = 3)
    
    add_curv_fit(x = myx, y = myy, 
                 withR2 = TRUE,
                 col ="gray80",
                 R2shiftX = -0.1, 
                 R2shiftY = -0.1,
                 lty=2)
    
    legend("topleft",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

}


################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
