library(foreach)
library(doMC)
library(tools)
library(DESeq2)

startTime <- Sys.time()

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 300, 7)

# Rscript cumul_gene_variance.R
# Rscript cumul_gene_variance.R <exprType>
# Rscript cumul_gene_variance.R fpkm
# Rscript cumul_gene_variance.R log2fpkm
# Rscript cumul_gene_variance.R NBvar
# Rscript cumul_gene_variance.R voom

# Rscript cumul_gene_variance.R log2fpkm 1000 raw



cat("> START: cumul_gene_variance.R\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) > 0)
exprType <- "log2fpkm"
if(!is.na(args[1]))
  exprType <- args[1]

nTopLast <- 1000

nTopColorDS <- 5

nValues_xvect <- 100 # # var_ref <- seq(from=0, to = overallMax, length.out = nValues_xvect) 

topCol <-  "darkorange2"
botCol <-  "dodgerblue4"
otherCol <- "gainsboro"


if(!is.na(args[2]))
  nTopLast <- as.numeric(args[2])

if(!is.na(args[3]))
  transfType <- args[3]

cat("... START with exprType =\t", exprType, "\n")
cat("... START with nTopLast =\t", nTopLast, "\n")
cat("... START with transfType =\t", transfType, "\n")

stopifnot(is.numeric(nTopLast))
stopifnot(exprType %in% c("fpkm", "log2fpkm", "NBvar", "voom"))
stopifnot(transfType %in% c("rescVar", "raw", "log10") | transfType %in% paste0(c("rescVar", "raw", "log10"), "_cropVar") )

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
stopifnot(file.exists(pipMainFolder))

outPipFold <- file.path(pipMainFolder, "OUTPUT_FOLDER")
stopifnot(file.exists(outPipFold))

aucFCCfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUCfcc_vs_AUCcoexprdist/all_auc_FCC.Rdata")
aucCoexprDistFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/AUCfcc_vs_AUCcoexprdist/all_auc_CoexprDist.Rdata")
stopifnot(file.exists(aucFCCfile))
stopifnot(file.exists(aucCoexprDistFile))

all_setting_files <- list.files(settingFolder, full.names=T)
all_setting_files <- all_setting_files[grep(".R$", all_setting_files)]
stopifnot(length(all_setting_files) > 0)

outFold <- file.path(paste0("CUMUL_GENE_VARIANCE"), paste0(toupper(exprType), "_", nTopLast), transfType)
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
# myHeight <- ifelse(plotType == "png", 600, 10)
# myWidth <- ifelse(plotType == "png", 600, 10)
myHeight <- ifelse(plotType == "png", 300, 5)
myWidth <- ifelse(plotType == "png", 300, 5)

# all_setting_files <- all_setting_files[1:3]

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
            meanMostVar = NA,
            meanLeastVar = NA,
            stringsAsFactors = FALSE
          ))
        }
          
        stopifnot(is.numeric(curr_exprDT[1,1]))
      }
    
    geneVar <- apply(curr_exprDT, 1,  var, na.rm=T)
    geneVar <- sort(geneVar, decreasing = TRUE)
      
    stopifnot(length(geneVar) >= nTopLast)
  
    keepTop <- c(1:nTopLast)
    mostVariant <- geneVar[keepTop]
    names(mostVariant) <- paste0(keepTop, "_", names(mostVariant))
    
    keepLast <- c((length(geneVar)-nTopLast+1):length(geneVar))
    leastVariant <- geneVar[keepLast]
    names(leastVariant) <- paste0(keepLast, "_", names(leastVariant))
    
    stopifnot(length(mostVariant) == nTopLast)
    stopifnot(length(leastVariant) == nTopLast)
    
    # meanMostVar <- mean(mostVariant)
    # meanLeastVar <- mean(leastVariant)
    
    data.frame(
      nTopLast = nTopLast,
      data_type = exprType,
      dataset = curr_ds,
      position = c(rep("mostVar", nTopLast), rep("leastVar", nTopLast)),
      gene = c(names(mostVariant), names(leastVariant)),
      var =  c(mostVariant, leastVariant),
      stringsAsFactors = FALSE
    )
    
  }
  rownames(all_ds_geneVarDT) <- NULL
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  save(all_ds_geneVarDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  all_ds_geneVarDT <- eval(parse(text = load(outFile)))
}



######################################################################################
###################################################################################### RETRIEVE FCC RANK
###################################################################################### 

all_ds <- unique(all_ds_geneVarDT$dataset)

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


#### DRAW THE CUMUL PLOT

myylab <- "% genes with variance >= threshold"

if(grepl("_cropVar", transfType)) {
  upLimit <- as.numeric(quantile(all_ds_geneVarDT$var, probs = 0.95))
  stopifnot(!is.na(upLimit))
  all_ds_geneVarDT$var[ all_ds_geneVarDT$var > upLimit ] <- upLimit
}

# if(transfType == "rescVar"){
if(grepl("rescVar", transfType)){
  all_ds_geneVarDT$var <- all_ds_geneVarDT$var/max(all_ds_geneVarDT$var)
  myxlab <- "Rescaled variance threshold"
# } else if(transfType == "log10"){
} else if(grepl( "log10", transfType)){  
  all_ds_geneVarDT$var <- log10(all_ds_geneVarDT$var + 1)
  myxlab <- "Variance threshold [log10(#+1)]"
# } else if(transfType == "raw"){
} else if(grepl("raw", transfType)){  
  myxlab <- "Variance threshold"
} else{
  stop("-- error -- \n")
}

overallMax <- max(all_ds_geneVarDT$var)
stopifnot(!is.na(overallMax))
var_ref <- seq(from=0, to = overallMax, length.out = nValues_xvect)

myxlim <- range(var_ref)

all_y <- list()

nDs <- length(unique(all_ds_geneVarDT$dataset))
mtextSubtit <- paste0("all DS (n=", nDs, ")")

# positionType="mostVar"
for(positionType in unique(all_ds_geneVarDT$position)) {
  
  # mySub <- paste0(positionType , " ", nTopLast, " - ", exprType)    
  
  mySub <- paste0("(", nTopLast, " genes with ", positionType, " - ",exprType, ")")
  
  sub_dt <- all_ds_geneVarDT[all_ds_geneVarDT$position == positionType,]
  all_y[[paste0(positionType)]] <- list()
  # curr_ds="TCGAcrc_msi_mss"
  for(curr_ds in unique(all_ds_geneVarDT$dataset)){
    curr_dt <- sub_dt[sub_dt$dataset == curr_ds,]
    stopifnot(nrow(curr_dt) > 0)
    
            # !!! WRONG should rescale for the overall variance !!!
            # # curr_dt$var_norm <- curr_dt$var/max(curr_dt$var)
            # if(rescVar){
            #   curr_dt$var_norm <- curr_dt$var/overallMax 
            # } else{
            #   curr_dt$var_norm <- curr_dt$var
            # }
    curr_dt$var_norm <- curr_dt$var
    stopifnot(curr_dt$var_norm >= min(var_ref) & curr_dt$var_norm <= max(var_ref))
    
    y_nGreaterVar <- sapply(var_ref, function(x) sum(curr_dt$var_norm >= x)/length(curr_dt$var_norm) )
    myTit <- paste0("% of genes with variance >= threshold") 
    
                # outFile <- file.path(outFold, paste0(positionType, "_", curr_ds, "_", transfType, ".", plotType))
                # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
                # plot(x = var_ref,
                #      y = round(y_nGreaterVar*100,2),
                #      type="l",
                #      xlab = myxlab,
                #      ylab = myylab,
                #      main = myTit)
                # mtext(side = 3, text = paste0(curr_ds))
                # mtext(side = 2, text = mySub, cex =0.7, line = 2)
                # 
                # foo <- dev.off()
                # cat(paste0("... written: ", outFile, "\n"))
    
    all_y[[paste0(positionType)]][[paste0(curr_ds)]] <- y_nGreaterVar
  }
  # myTit <- paste0("Fraction of genes with variance >= threshold\nall ds") 
  myTit <- paste0("% of genes with variance >= threshold") 
  
  
                      # outFile <- file.path(outFold, paste0(positionType, "_allDS_", transfType, ".", plotType))
                      # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
                      # plot(NULL,
                      #      # xlim = c(0,1),
                      #      xlim = myxlim,
                      #      # ylim = c(0,1),
                      #      ylim = c(0,100),
                      #      xlab = myxlab,
                      #      ylab = myylab,
                      #      main = myTit)
                      # mtext(side = 3, text = mtextSubtit)   
                      # mtext(side = 2, text = mySub, cex =0.7, line = 2)
                      # 
                      # foo <- lapply(all_y[[paste0(positionType)]], function(x) lines(x=var_ref, y=round(x*100,2),
                      #                                                         col = "black"))
                      # foo <- dev.off()
                      # cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(positionType, "_allDS_botTopColors_", transfType, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(NULL,
       # xlim = c(0,1),
       xlim = myxlim,
       # ylim = c(0,1),
       ylim = c(0,100),
       xlab = myxlab,
       ylab = myylab,
       main = myTit)
  mtext(side = 3, text = mtextSubtit)   
  mtext(side = 2, text = mySub, cex =0.7, line = 2)
  
  # foo <- lapply(all_y[[paste0(positionType)]], function(x) lines(x=var_ref, y=round(x*100,2),
  #                                                                col = "black"))
  # 
  foo <- lapply(1:length(all_y[[paste0(positionType)]]), function(i) {
    dsname <- names(all_y[[paste0(positionType)]])[i]
    x <- unlist(all_y[[paste0(positionType)]][i])
    lineCol <- ifelse(dsname %in% topDS, topCol,
                      ifelse(dsname %in% botDS, botCol, otherCol))
    lines(x=var_ref, y=round(x*100,2),
          col = lineCol)
  })

  foo <- lapply(1:length(all_y[[paste0(positionType)]]), function(i) {
    dsname <- names(all_y[[paste0(positionType)]])[i]
    x <- unlist(all_y[[paste0(positionType)]][i])
    lineCol <- ifelse(dsname %in% topDS, topCol,
                      ifelse(dsname %in% botDS, botCol, otherCol))
    if(dsname %in% c(topDS, botDS)) {
      lines(x=var_ref, y=round(x*100,2),
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


mySub <- paste0(nTopLast, " - ", exprType)    

outFile <- file.path(outFold, paste0("mostVar_leastVar_allDS_", transfType, ".", plotType))
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

foo <- lapply(all_y[[paste0(unique(all_ds_geneVarDT$position)[1])]], function(x) lines(x=var_ref, y=round(x*100,2),
                                                        col = "blue"))
foo <- lapply(all_y[[paste0(unique(all_ds_geneVarDT$position)[2])]], function(x) lines(x=var_ref, y=round(x*100,2),
                                                           col = "red"))
legend("topright", lty=1, legend = c(paste0(unique(all_ds_geneVarDT$position)[1]), paste0(unique(all_ds_geneVarDT$position)[2])),
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
