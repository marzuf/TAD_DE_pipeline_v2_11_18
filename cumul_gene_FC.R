library(foreach)
library(doMC)

startTime <- Sys.time()

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 300, 7)

# Rscript cumul_gene_FC.R [<nTopLast>] [<transfType>]
# Rscript cumul_gene_FC.R 1000 rescFC
# Rscript cumul_gene_FC.R 1000 raw_cropFC

cat("> START: cumul_gene_FC.R\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18")

args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) > 0)


nTopLast <- 1000

nTopColorDS <- 5

nValues_xvect <- 100 # # fc_ref <- seq(from=0, to = overallMax, length.out = 100) 

topCol <-  "darkorange2"
botCol <-  "dodgerblue4"
otherCol <- "gainsboro"

transfType <- "rescFC"

if(!is.na(args[1])) {
  nTopLast <- as.numeric(args[1])
  stopifnot(is.numeric(nTopLast))
}
if(!is.na(args[2])) {
  transfType <- args[2]
}
stopifnot(transfType %in% c("rescFC", "raw", "log10") | transfType %in% paste0(c("rescFC", "raw", "log10"), "_cropFC") )

cat("... START with nTopLast =\t", nTopLast, "\n")
cat("... START with transfType =\t", transfType, "\n")

buildTable <- T

# registerDoMC(ifelse(SSHFS, 2, 20))

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))

settingFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

pipMainFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")
stopifnot(file.exists(pipMainFolder))

outPipFold <- file.path(pipMainFolder, "OUTPUT_FOLDER")
stopifnot(file.exists(outPipFold))

all_setting_files <- list.files(settingFolder, full.names=T)
all_setting_files <- all_setting_files[grep(".R$", all_setting_files)]
stopifnot(length(all_setting_files) > 0)

outFold <- file.path(paste0("CUMUL_GENE_FC"), paste0(nTopLast), transfType)
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
# myHeight <- ifelse(plotType == "png", 600, 10)
# myWidth <- ifelse(plotType == "png", 600, 10)
myHeight <- ifelse(plotType == "png", 300, 5)
myWidth <- ifelse(plotType == "png", 300, 5)

ds_file="/media/electron//mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput/run_settings_TCGAcrc_msi_mss.R"
ds_file="/media/electron//mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput/run_settings_GSE102073_stic_nostic.R"
if(buildTable) {
  all_ds_geneFC_DT <- foreach(ds_file = all_setting_files, .combine="rbind") %dopar% {
  
    stopifnot(file.exists(ds_file))
    cat("... source settingFile", basename(ds_file), "\n")
    source(ds_file)
  
    curr_ds <- basename(pipOutFold)
    cat("... START", curr_ds, "\n")
  
    ds_pipFolder <- file.path(pipMainFolder, pipOutFold)
    cat(ds_pipFolder,"\n")
    stopifnot(file.exists(ds_pipFolder))
    
    
    # load the fold change values
    cat("... load limma DE results\n")
    limmaDT <- eval(parse(text = load(
      file.path(ds_pipFolder, "1_runGeneDE", paste0("DE_topTable.Rdata")))))
        
    cat("... load geneList\n")
    geneList <- eval(parse(text = load(
      file.path(ds_pipFolder, "0_prepGeneData", "pipeline_geneList.Rdata")
    )))
    
    
    # retain only the genes used in the pipeline
    stopifnot(rownames(limmaDT) == limmaDT$genes)
    stopifnot(names(geneList) %in% rownames(limmaDT))
    
    sub_limmaDT <- limmaDT[names(geneList),]
    stopifnot(nrow(sub_limmaDT) == length(geneList))
    stopifnot(rownames(sub_limmaDT) == names(geneList))
    stopifnot(rownames(sub_limmaDT) == sub_limmaDT$genes)
    stopifnot(sub_limmaDT$genes == names(geneList))
    
    genes_absFC <- abs(setNames(sub_limmaDT$logFC, sub_limmaDT$genes))
    genes_absFC_sort <- sort(genes_absFC, decreasing = TRUE)
    
    stopifnot(length(genes_absFC_sort) >= nTopLast)
  
    keepTop <- c(1:nTopLast)
    mostFC <- genes_absFC_sort[keepTop]
    names(mostFC) <- paste0(keepTop, "_", names(mostFC))
    
    keepLast <- c((length(genes_absFC_sort)-nTopLast+1):length(genes_absFC_sort))
    leastFC <- genes_absFC_sort[keepLast]
    names(leastFC) <- paste0(keepLast, "_", names(leastFC))
    
    stopifnot(length(mostFC) == nTopLast)
    stopifnot(length(leastFC) == nTopLast)
    
    data.frame(
      nTopLast = nTopLast,
      dataset = curr_ds,
      position = c(rep("mostFC", nTopLast), rep("leastFC", nTopLast)),
      gene = c(names(mostFC), names(leastFC)),
      absFC =  c(mostFC, leastFC),
      stringsAsFactors = FALSE
    )
    
  }
  rownames(all_ds_geneFC_DT) <- NULL
  outFile <- file.path(outFold, "all_ds_geneFC_DT.Rdata")
  save(all_ds_geneFC_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_geneFC_DT.Rdata")
  all_ds_geneFC_DT <- eval(parse(text = load(outFile)))
}


#### DRAW THE CUMUL PLOT

range(all_ds_geneFC_DT$absFC)


if(grepl("_cropFC", transfType)) {
  upLimit <- as.numeric(quantile(all_ds_geneFC_DT$absFC, probs = 0.95))
  stopifnot(!is.na(upLimit))
  all_ds_geneFC_DT$absFC[ all_ds_geneFC_DT$absFC > upLimit ] <- upLimit
}

# if(transfType == "rescFC") {
if( grepl("rescFC", transfType)) {
  myxlab <- "Rescaled FC threshold" 
  all_ds_geneFC_DT$absFC <- all_ds_geneFC_DT$absFC/max(all_ds_geneFC_DT$absFC)
# else if(transfType == "rescFC") {
} else if( grepl("log10", transfType)) {
  myxlab <- "FC threshold [log10(#+1)]"
  all_ds_geneFC_DT$absFC <- log10(all_ds_geneFC_DT$absFC+1)
# } else if(transfType == "raw") {
} else if(grepl("raw", transfType)) {
  myxlab <- "FC threshold"
} else{
  stop("-- error --\n")
}

overallMax <- max(all_ds_geneFC_DT$absFC)
cat(overallMax,"\n")

stopifnot(!is.na(overallMax))
fc_ref <- seq(from=0, to = overallMax, length.out = nValues_xvect) 

myylab_init <- "% genes with FC >= threshold"

myxlim <- range(fc_ref)
all_y <- list()

range(all_ds_geneFC_DT$absFC)


######################################################################################
###################################################################################### RETRIEVE FCC RANK
###################################################################################### 

all_ds <- unique(all_ds_geneFC_DT$dataset)

all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(outPipFold, curr_ds, "170_score_auc_pval_withShuffle")
  stopifnot(file.exists(step17_fold))
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  stopifnot(file.exists(aucFCC_file))
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  aucFCC
}


names(all_ds_aucFCC) <- all_ds
all_ds_aucFCC <- sort(all_ds_aucFCC, decreasing=TRUE)
topDS <- names(all_ds_aucFCC[1:nTopColorDS])
botDS <- names(rev(all_ds_aucFCC)[1:nTopColorDS])

######################################################################################
###################################################################################### PLOT SEPARATELY topLeast AND topMost; 1 DS + ALL DS
###################################################################################### 

cat(range(all_ds_geneFC_DT$absFC),"\n")


# positionType="mostFC"
for(positionType in unique(all_ds_geneFC_DT$position)) {
  # mySub <- paste0(positionType , " ", nTopLast)    
  mySub <- paste0("(", nTopLast, " genes with ", positionType, ")")
  
  # myylab <- paste0(myylab_init, "\n", mySub)
  myylab <- paste0(myylab_init)
  
  sub_dt <- all_ds_geneFC_DT[all_ds_geneFC_DT$position == positionType,]
  all_y[[paste0(positionType)]] <- list()
  
  nDs <- length(unique(all_ds_geneFC_DT$dataset))
  mtextSubtit <- paste0("all DS (n=", nDs, ")")
  
  # curr_ds="TCGAcrc_msi_mss"
  for(curr_ds in unique(all_ds_geneFC_DT$dataset)){
    curr_dt <- sub_dt[sub_dt$dataset == curr_ds,]
    stopifnot(nrow(curr_dt) > 0)
    ### !!! WRONG !!! SHOULD DIVIDE BY THE OVERALL MAX SO THAT I CAN COMPARE ACROSS DATASETS
    # curr_dt$fc_norm <- curr_dt$absFC/max(curr_dt$absFC)
    # if(transfType == "rescFC") {
    #   curr_dt$fc_norm <- curr_dt$absFC/overallMax
    # } else{
    #   
    # }
    curr_dt$fc_norm <- curr_dt$absFC
    
    stopifnot(curr_dt$fc_norm >= min(fc_ref) & curr_dt$fc_norm <= max(fc_ref))
    y_nGreaterFC <- sapply(fc_ref, function(x) sum(curr_dt$fc_norm >= x)/length(curr_dt$fc_norm) )
    
    myTit <- paste0("% of genes with FC >= threshold\n", curr_ds) 
    
                    # outFile <- file.path(outFold, paste0(positionType, "_", curr_ds, "_", transfType, ".", plotType))
                    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
                    # plot(x = fc_ref,
                    #      y = round(y_nGreaterFC*100,2),
                    #      type="l",
                    #      xlab = myxlab,
                    #      ylab = myylab,
                    #      # sub = mySub,
                    #      main = myTit
                    #      )
                    # mtext(side = 3, text = curr_ds)
                    # mtext(side = 2, text = mySub, cex =0.7)
                    # 
                    # foo <- dev.off()
                    # cat(paste0("... written: ", outFile, "\n"))
    
    all_y[[paste0(positionType)]][[paste0(curr_ds)]] <- y_nGreaterFC
  }
  myTit <- paste0("% of genes with FC >= threshold") 
  
                                      # outFile <- file.path(outFold, paste0(positionType, "_allDS_",  transfType, ".", plotType))
                                      # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
                                      # plot(NULL,
                                      #      # xlim = c(0,1),
                                      #      xlim = myxlim,
                                      #      # ylim = c(0,1),
                                      #      ylim = c(0,100),
                                      #      xlab = myxlab,
                                      #      ylab = myylab,
                                      #      # sub = mySub,
                                      #      main = myTit)
                                      # mtext(side = 3, text = mtextSubtit)   
                                      # mtext(side = 2, text = mySub, cex =0.7)
                                      # 
                                      # foo <- lapply(all_y[[paste0(positionType)]], function(x) lines(x=fc_ref, y=round(x*100,2),
                                      #                                                         col = "black"))
                                      # foo <- dev.off()
                                      # cat(paste0("... written: ", outFile, "\n"))
                                      # 
  
      # COLOR THE TOP AND BOTTOM DATASETS CURVES
      outFile <- file.path(outFold, paste0(positionType, "_allDS_botTopColors_",transfType,  ".", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot(NULL,
           # xlim = c(0,1),
           xlim = myxlim,
           # ylim = c(0,1),
           ylim = c(0,100),
           xlab = myxlab,
           ylab = myylab,
           # sub = mySub,
           main = myTit
           )
      mtext(side = 3, text = mtextSubtit)
      mtext(side = 2, text = mySub, cex =0.7, line = 2)

      # foo <- lapply(all_y[[paste0(positionType)]], function(x) {
        foo <- lapply(1:length(all_y[[paste0(positionType)]]), function(i) {
          dsname <- names(all_y[[paste0(positionType)]])[i]
          x <- unlist(all_y[[paste0(positionType)]][i])
          lineCol <- ifelse(dsname %in% topDS, topCol,
                              ifelse(dsname %in% botDS, botCol, otherCol))
          lines(x=fc_ref, y=round(x*100,2),
              col = lineCol)
      })

        foo <- lapply(1:length(all_y[[paste0(positionType)]]), function(i) {
          dsname <- names(all_y[[paste0(positionType)]])[i]
          x <- unlist(all_y[[paste0(positionType)]][i])
          lineCol <- ifelse(dsname %in% topDS, topCol,
                              ifelse(dsname %in% botDS, botCol, otherCol))
          if(dsname %in% c(topDS, botDS)) {
            lines(x=fc_ref, y=round(x*100,2),
              col = lineCol)
          }
      })

      legend("topright", lty=1, legend = c(paste0("top ", nTopColorDS, " FCC DS"), paste0("last ", nTopColorDS, " FCC DS")),
                                           col=c(topCol, botCol), bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
}

######################################################################################
###################################################################################### 1 PLOTH ALL DS +  topLeast AND topMost
###################################################################################### 

mySub <- paste0(nTopLast)    

outFile <- file.path(outFold, paste0("mostFC_leastFC_allDS_",transfType, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     # xlim = c(0,1),
     xlim = myxlim,
     # ylim = c(0,1),
     ylim = c(0,100),
     xlab = myxlab,
     ylab = myylab,
     main = myTit)
mtext(side = 3, text = mySub)   
foo <- lapply(all_y[[paste0(unique(all_ds_geneFC_DT$position)[1])]], function(x) lines(x=fc_ref, y=round(x*100,2),
                                                        col = "blue"))
foo <- lapply(all_y[[paste0(unique(all_ds_geneFC_DT$position)[2])]], function(x) lines(x=fc_ref, y=round(x*100,2),
                                                           col = "red"))
legend("topright", lty=1, legend = c(paste0(unique(all_ds_geneFC_DT$position)[1]), paste0(unique(all_ds_geneFC_DT$position)[2])),
       col = c("blue", "red"),
       bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################################################################
outFile <- file.path(outFold, paste0("all_y_", transfType, ".Rdata"))
save(all_y, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
