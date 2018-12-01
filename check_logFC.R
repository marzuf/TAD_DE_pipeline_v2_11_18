startTime <- Sys.time()
cat(paste0("> Rscript check_coexpr_dist.R\n"))

# Rscript check_coexpr_dist.R

source("analysis_utils.R")

library(foreach)
library(doMC)

# source("coreg_utils.R")

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
pipOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")
all_datasets <- list.files(pipOutFold)
# head(all_datasets)
stopifnot(length(all_datasets) >  0 )
cat(paste0("... found ", length(all_datasets), " datasets\n"))

sizeLimit <- 500 * 1000

geneVarianceType <- "GENE_VARIANCE_LOG2FPKM"


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)


outFold <- file.path("CHECK_FC_VAR")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "check_coexpr_logFile.txt")
system(paste0("rm -f ", logFile))

corr_method <- "pearson"

# all_datasets = all_datasets[1]
# all_datasets = all_datasets[1:5]

load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18", "CUMUL_GENE_FC/1000/raw/all_ds_geneFC_DT.Rdata"))
all_ds_geneFC_DT <- all_ds_geneFC_DT[order(all_ds_geneFC_DT$absFC, decreasing = T),]
head(all_ds_geneFC_DT)

all_ds_geneFC_DT_mean <- aggregate(absFC ~ dataset, FUN=mean, data = all_ds_geneFC_DT)
dataset_meanFC <- setNames(all_ds_geneFC_DT_mean$absFC, all_ds_geneFC_DT_mean$dataset)

all_ds_geneFC_DT_max <- aggregate(absFC ~ dataset, FUN=max, data = all_ds_geneFC_DT)
dataset_maxFC <- setNames(all_ds_geneFC_DT_max$absFC, all_ds_geneFC_DT_max$dataset)



load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18", "CUMUL_GENE_VARIANCE/LOG2FPKM_1000/raw/all_ds_geneVarDT.Rdata"))
all_ds_geneVarDT <- all_ds_geneVarDT[order(all_ds_geneVarDT$var, decreasing = T),]
head(all_ds_geneVarDT)
all_ds_geneVarDT_mean <- aggregate(var ~ dataset, FUN=mean, data = all_ds_geneVarDT)
dataset_meanVar <- setNames(all_ds_geneVarDT_mean$var, all_ds_geneVarDT_mean$dataset)

all_ds_geneVarDT_max <- aggregate(var ~ dataset, FUN=max, data = all_ds_geneVarDT)
dataset_maxVar <- setNames(all_ds_geneVarDT_max$var, all_ds_geneVarDT_max$dataset)



if(buildTable) {
  
  datasets_coexpr_summary_DT <- foreach(curr_dataset = all_datasets, .combine='rbind') %dopar% {
    
    # txt <- paste0("> START ", curr_dataset, "\n")
    # printAndLog(txt, logFile)
    
    cat("... load allData\n")
    # data_dt <- eval(parse(text = load(file.path(setDir, 
    #                                             "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUC_COEXPRDIST_SORTNODUP",curr_dataset, "allData_dt.Rdata"))))
    # head(data_dt)
    # 
    # cat("... load qqnorm\n")
    # qqnorm_dt <- eval(parse(text = load(file.path(setDir,
    #                                               "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER", curr_dataset, "0_prepGeneData/rna_qqnorm_rnaseqDT.Rdata"))))
    # dim(qqnorm_dt)
    # qqnorm_dt[1:5,1:5]
    
    cat("... load rnaseq\n")
    # rnaseq_dt <- eval(parse(text = load(file.path(setDir,
    #                                               "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER", curr_dataset, "0_prepGeneData/rna_rnaseqDT.Rdata"))))
    # dim(rnaseq_dt)
    # rnaseq_dt[1:5,1:5]
    
    # cat("... load coexpr\n")
    # coexpr_dt <- eval(parse(text = load(file.path(setDir,
    #                                               "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_COEXPR_SORTNODUP", paste0(curr_dataset, "_", corr_method), "coexprDT.Rdata"
    # ))))
    # dim(coexpr_dt)
    # head(coexpr_dt)
    
    # cat("... load pipeline_geneList\n")
    # pipeline_geneList <- eval(parse(text = load(file.path(setDir,
    #                                                       "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER", curr_dataset, "0_prepGeneData/pipeline_geneList.Rdata"))))
    # head(pipeline_geneList)
    # 
    # g1 <- as.character(data_dt$gene1[1])
    # g2 <- as.character(data_dt$gene2[1])
    # 
    # g1_match <- names(pipeline_geneList)[as.character(pipeline_geneList) == g1]
    # g2_match <- names(pipeline_geneList)[as.character(pipeline_geneList) == g2]
    
    # expr1 <- rnaseq_dt[g1_match,]
    # expr2 <- rnaseq_dt[g2_match,]
    
    # qqexpr1 <- qqnorm_dt[g1_match,]
    # qqexpr2 <- qqnorm_dt[g2_match,]
    # 
    # cor_expr <- cor(as.numeric(expr1), as.numeric(expr2), method = corr_method)
    # cor_expr_spearman <- cor(as.numeric(expr1), as.numeric(expr2), method = "spearman")
    # cor_exprQQ <- cor(as.numeric(qqexpr1), as.numeric(qqexpr2), method = corr_method)
    # cor_exprQQ_spearman <- cor(as.numeric(qqexpr1), as.numeric(qqexpr2), method = "spearman")
    
    # calc_coexpr <- coexpr_dt$coexpr[as.character(coexpr_dt$gene1) == g1 & as.character(coexpr_dt$gene2) == g2]
    
    # stopifnot(round(cor_exprQQ, 4) == round(calc_coexpr,4))
    
    mytit <- curr_dataset
    mysubtit <- paste0("sizeLimit = ", sizeLimit/1000, "kb")
    
    # outFile <- file.path(outFold, paste0(curr_dataset, "_density_limit_", sizeLimit/1000,"kb.", plotType))
    # do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    # plot(density(data_dt$coexpr[data_dt$dist <= sizeLimit]),
    #      main = mytit)
    # mtext(text = mysubtit, side=3)
    # cat(paste0("... written: ", outFile, "\n"))
    # foo <- dev.off()
    # 
    settingF <- file.path(setDir,
                          "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput",
                          paste0("run_settings_", curr_dataset, ".R"))
    stopifnot(file.exists(settingF))
    source(settingF)
    sample1_file <- file.path(setDir, sample1_file)
    stopifnot(file.exists(sample1_file))
    s1 <- eval(parse(text = load(sample1_file)))
    
    sample2_file <- file.path(setDir, sample2_file)
    stopifnot(file.exists(sample2_file))
    s2 <- eval(parse(text = load(sample2_file)))
    
    # txt <- paste0("... cor_expr =\t", round(cor_expr, 4), "\n")
    # printAndLog(txt, logFile)
    # 
    # txt <- paste0("... cor_exprQQ =\t", round(cor_exprQQ, 4), "\n")
    # printAndLog(txt, logFile)
    # 
    # txt <- paste0("... calc_coexpr =\t", round(calc_coexpr, 4), "\n")
    # printAndLog(txt, logFile)
    # 
    # txt <- paste0("... cor_expr_spearman =\t", round(cor_expr_spearman, 4), "\n")
    # printAndLog(txt, logFile)
    # 
    # txt <- paste0("... cor_exprQQ_spearman =\t", round(cor_exprQQ_spearman, 4), "\n")
    # printAndLog(txt, logFile)
    # 
    # 
    # sink(file = logFile, append =T)
    # print(summary(data_dt$coexpr[data_dt$dist < sizeLimit]))
    # cat("\n")
    # sink()
    
    # summary_coexpr <- summary(data_dt$coexpr[data_dt$dist <= sizeLimit])
    
    data.frame(
      
      dataset=curr_dataset,
      nSamp1 = length(s1),
      nSamp2 = length(s2),
      geneID_format = geneID_format,
      inputDataType = inputDataType,
      # nGenePairs = length(data_dt$coexpr[data_dt$dist <= sizeLimit]),
      # coexpr_mean = as.numeric(summary_coexpr["Mean"]),
      # coexpr_median = as.numeric(summary_coexpr["Median"]),
      # coexpr_firstQ = as.numeric(summary_coexpr["1st Qu."]),
      # coexpr_thirdQ = as.numeric(summary_coexpr["3rd Qu."]),
      meanVar = dataset_meanVar[curr_dataset],
      meanFC = dataset_meanFC[curr_dataset],
      maxVar = dataset_maxVar[curr_dataset],
      maxFC = dataset_maxFC[curr_dataset],
      stringsAsFactors = FALSE
      
    )
    
  }
  
  outFile <- file.path(outFold, "datasets_coexpr_summary_DT.Rdata")
  save(datasets_coexpr_summary_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFold, "datasets_coexpr_summary_DT.Rdata")
  stopifnot(file.exists(outFile))
  datasets_coexpr_summary_DT <- eval(parse(text = load(outFile)))
}

datasets_coexpr_summary_DT <- datasets_coexpr_summary_DT[order(datasets_coexpr_summary_DT$meanFC, decreasing=TRUE),]

datasets_coexpr_summary_DT_out <- datasets_coexpr_summary_DT
# datasets_coexpr_summary_DT_out$coexpr_mean <- round(datasets_coexpr_summary_DT_out$coexpr_mean, 4)
# datasets_coexpr_summary_DT_out$coexpr_median <- round(datasets_coexpr_summary_DT_out$coexpr_median, 4)
# datasets_coexpr_summary_DT_out$coexpr_firstQ <- round(datasets_coexpr_summary_DT_out$coexpr_firstQ, 4)
# datasets_coexpr_summary_DT_out$coexpr_thirdQ <- round(datasets_coexpr_summary_DT_out$coexpr_thirdQ, 4)

datasets_coexpr_summary_DT_out$meanVar <- round(datasets_coexpr_summary_DT_out$meanVar, 4)
datasets_coexpr_summary_DT_out$meanFC <- round(datasets_coexpr_summary_DT_out$meanFC, 4)

outFile <- file.path(outFold, "datasets_coexpr_summary_DT.txt")
write.table(datasets_coexpr_summary_DT_out, file = outFile, quote=F, sep="\t", col.names=T, row.names=F, append=F)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(nrow(datasets_coexpr_summary_DT) == length(all_datasets))

datasets_coexpr_summary_DT$nAllSamp <- datasets_coexpr_summary_DT$nSamp1 + datasets_coexpr_summary_DT$nSamp2 

all_ref_var <- c("meanVar", "meanFC", "maxVar", "maxFC")
all_curr_var <- c("nAllSamp")

# load the variance
varianceFile <- file.path(setDir, 
                          "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                          geneVarianceType,
                          "all_ds_geneVarDT.Rdata")
stopifnot(file.exists(varianceFile))

varianceDT <- eval(parse(text = load(varianceFile)))


stopifnot(sort(datasets_coexpr_summary_DT$dataset) == sort(varianceDT$dataset))
datasets_coexpr_summary_DT <- merge(datasets_coexpr_summary_DT, varianceDT, by="dataset")
stopifnot(sort(datasets_coexpr_summary_DT$dataset) == sort(varianceDT$dataset))

for(ref_var in all_ref_var) {
  
  for(curr_var in all_curr_var) {
    
    myx <- datasets_coexpr_summary_DT[, curr_var]
    myy <- datasets_coexpr_summary_DT[,ref_var]
    
    # outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var,  "_limit_", sizeLimit/1000,"kb.", plotType))
    outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var, ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(x = myx,
         y = myy,
         pch=16, cex=0.7,
         xlab = curr_var,
         ylab = ref_var,
        main = paste0(ref_var, " vs. ", curr_var))    
    mtext(paste0("(n = ", nrow(datasets_coexpr_summary_DT),")"), side = 3)
    add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
    addCorr(x=myx, y=myy, legPos="topright", corMet="pearson", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var,  "_limit_", sizeLimit/1000,"kb_withLabs.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(x = myx,
         y = myy,
         pch=16, cex=0.7,
         xlab = curr_var,
         ylab = ref_var,
        main = paste0(ref_var, " vs. ", curr_var))    
    text(x = myx,
         y = myy,
         labels = datasets_coexpr_summary_DT$dataset,
         cex=0.7
         )
    mtext(paste0("(n = ", nrow(datasets_coexpr_summary_DT),")"), side = 3)
    add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
    addCorr(x=myx, y=myy, legPos="topright", corMet="pearson", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}




curr_var = "meanFC"
ref_var = "meanVar"


myx <- datasets_coexpr_summary_DT[, curr_var]
myy <- datasets_coexpr_summary_DT[,ref_var]

outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var,  "_limit_", sizeLimit/1000,"kb.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = myx,
     y = myy,
     pch=16, cex=0.7,
     xlab = curr_var,
     ylab = ref_var,
     main = paste0(ref_var, " vs. ", curr_var))    
mtext(paste0("(n = ", nrow(datasets_coexpr_summary_DT),")"), side = 3)
add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
addCorr(x=myx, y=myy, legPos="topright", corMet="pearson", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var,  "_limit_", sizeLimit/1000,"kb_withLabs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = myx,
     y = myy,
     pch=16, cex=0.7,
     xlab = curr_var,
     ylab = ref_var,
     main = paste0(ref_var, " vs. ", curr_var))    
text(x = myx,
     y = myy,
     labels = datasets_coexpr_summary_DT$dataset,
     cex=0.7
)
mtext(paste0("(n = ", nrow(datasets_coexpr_summary_DT),")"), side = 3)
add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
addCorr(x=myx, y=myy, legPos="topright", corMet="pearson", bty="n")
foo <- dev.off()


curr_var = "maxFC"
ref_var = "maxVar"


myx <- datasets_coexpr_summary_DT[, curr_var]
myy <- datasets_coexpr_summary_DT[,ref_var]

outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var,  "_limit_", sizeLimit/1000,"kb.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = myx,
     y = myy,
     pch=16, cex=0.7,
     xlab = curr_var,
     ylab = ref_var,
     main = paste0(ref_var, " vs. ", curr_var))    
mtext(paste0("(n = ", nrow(datasets_coexpr_summary_DT),")"), side = 3)
add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
addCorr(x=myx, y=myy, legPos="topright", corMet="pearson", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("all_datasets_", ref_var, "_vs_", curr_var,  "_limit_", sizeLimit/1000,"kb_withLabs.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = myx,
     y = myy,
     pch=16, cex=0.7,
     xlab = curr_var,
     ylab = ref_var,
     main = paste0(ref_var, " vs. ", curr_var))    
text(x = myx,
     y = myy,
     labels = datasets_coexpr_summary_DT$dataset,
     cex=0.7
)
mtext(paste0("(n = ", nrow(datasets_coexpr_summary_DT),")"), side = 3)
add_curv_fit(x=myx, y=myy, withR2 = FALSE,lty=2)
addCorr(x=myx, y=myy, legPos="topright", corMet="pearson", bty="n")
foo <- dev.off()


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
