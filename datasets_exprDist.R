# Rscript datasets_exprDist.R

cat("> START: gene_variance.R\n")

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

cancerDS <- score_DT$dataset[score_DT$process_short == "cancer"]
noTCGA_cancerDS <- cancerDS[!grepl("^TCGA", cancerDS)]
TCGA_cancerDS <- cancerDS[grepl("TCGA", cancerDS)]
no_cancerDS <- score_DT$dataset[score_DT$process_short != "cancer"]

cat("... found # noTCGA_cancerDS\t=\t", length(noTCGA_cancerDS) , "\n" )
cat("... found # TCGA_cancerDS\t=\t", length(TCGA_cancerDS) , "\n" )
cat("... found # no_cancerDS\t=\t", length(no_cancerDS) , "\n" )

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

buildTable <- TRUE

# registerDoMC(ifelse(SSHFS, 2, 20))

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))

settingFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

pipMainFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

dsFold <- file.path(pipMainFolder, "OUTPUT_FOLDER")

all_setting_files <- list.files(settingFolder, full.names=T)
all_setting_files <- all_setting_files[grep(".R$", all_setting_files)]
stopifnot(length(all_setting_files) > 0)

outFold <- file.path(paste0("DATASETS_EXPRDIST"), toupper(exprType))
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

# all_setting_files <- all_setting_files[1:5]

if(buildTable) {
  all_ds_fpkmCounts <- foreach(ds_file = all_setting_files) %dopar% {
    
    if(!file.exists(ds_file)) return(NULL)
    stopifnot(file.exists(ds_file))
    cat("... source settingFile", basename(ds_file), "\n")
    source(ds_file)
    
    curr_ds <- basename(pipOutFold)
    cat("... START", curr_ds, "\n")
    
    ds_pipFolder <- file.path(pipMainFolder, pipOutFold)
    if(!file.exists(ds_pipFolder)) return(NULL)
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
    
    bpCounts <- boxplot(as.numeric(data.matrix(curr_exprDT)), 
                        plot=FALSE)
    
    dataType <- inputDataType
    
    list(bpCounts = bpCounts, dataType = dataType)
  }
  dataset_names <- gsub("run_settings_(.+)\\.R", "\\1", basename(all_setting_files))
  names(all_ds_fpkmCounts) <- all_setting_files
  # names(all_ds_fpkmCounts) <- all_setting_files
  outFile <- file.path(outFold, "all_ds_fpkmCounts.Rdata")
  save(all_ds_fpkmCounts, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_fpkmCounts.Rdata")
  #outFile <- "DATASETS_EXPRDIST/LOG2FPKM/all_ds_fpkmCounts.Rdata"
  all_ds_fpkmCounts <- eval(parse(text = load(outFile)))
}
stop("--ok\n")
# all_ds_fpkmCounts$color <- ifelse(all_ds_fpkmCounts$dataset %in% noTCGA_cancerDS, "green",
#                                   ifelse(all_ds_fpkmCounts$dataset %in% TCGA_cancerDS, "red",
#                                          ifelse(all_ds_fpkmCounts$dataset %in% no_cancerDS, "black", NA)))
# stopifnot(!is.na(all_ds_fpkmCounts$color))
# 
# 
# all_ds_fpkmCounts$dsType <- ifelse(all_ds_fpkmCounts$dataset %in% noTCGA_cancerDS, "cancer_notTCGA",
#                                    ifelse(all_ds_fpkmCounts$dataset %in% TCGA_cancerDS, "TCGA",
#                                           ifelse(all_ds_fpkmCounts$dataset %in% no_cancerDS, "not_cancer", NA)))
# 
# stopifnot(!is.na(all_ds_fpkmCounts$dsType))

boxplot_y_max <- max(unlist(lapply(all_ds_fpkmCounts, function(x) {
  bxpVar <- x[["bpCounts"]]
  max(bxpVar$out, na.rm=T)
})))

boxplot_y_min <- min(unlist(lapply(all_ds_fpkmCounts, function(x) {
  bxpVar <- x[["bpCounts"]][["out"]]
  min(bxpVar$out, na.rm=T)
})))


plot(NULL, xlab="", ylab="", main="",
     ylim = c(boxplot_y_min, boxplot_y_max),
     xlim = c(0, length(all_ds_fpkmCounts))
)

i=1
for(i in seq_len(nrow(all_ds_fpkmCounts))) {
  
  curr_ds <- names(all_ds_fpkmCounts)[i]
  
  curr_col  <- ifelse(curr_ds %in% noTCGA_cancerDS, "green",
                                        ifelse(curr_ds %in% TCGA_cancerDS, "red",
                                               ifelse(curr_ds %in% no_cancerDS, "black", NA)))
  stopifnot(!is.na(curr_col))
  
  
  # retrieve the boxplot
  boxplotInfo <- all_ds_fpkmCounts[[i]][["bpCounts"]]
  dsType <- all_ds_fpkmCounts[[i]][["dataType"]]
  
  bxp(boxplotInfo, add = TRUE, col=curr_col)
  
  
}

################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))