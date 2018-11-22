startTime <- Sys.time()
cat(paste0("> Rscript topTADs_density_withRandom.R\n"))

# Rscript topTADs_density_withRandom.R

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(reshape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

source("analysis_utils.R")

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "png"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 500, 10)
myHeightGG <- 7
myWidthGG <- 10

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18")

nTopTADs <- 20
# nTopTADs <- "all"
nToPlot <- 5
nRandom <- 1000

stopifnot(nTopTADs > 0)
stopifnot(nToPlot > 0)
stopifnot(nRandom > 0)

outFold <- file.path("TOP_TADs_DENSITY_WITH_RANDOM",  paste0(nTopTADs, "_rand", nRandom))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "top_tads_density.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"

all_ds <- list.files(dsFold)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

mySub <- paste0("(nTopTADs = ", nTopTADs , ")") 

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... nTopTADs:\t", nTopTADs, "\n")
printAndLog(txt, logFile)
txt <- paste0("... nToPlot:\t", nToPlot, "\n")
printAndLog(txt, logFile)


if(nTopTADs!="all") stopifnot(nTopTADs > 0)


# all_ds <- all_ds[1:2]

#all_ds <- all_ds[1:3]
curr_ds="TCGAcrc_msi_mss"
if(buildTable){
  all_ds_DT <- foreach(curr_ds = all_ds, .combine='rbind') %do% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    printAndLog(txt, logFile)
    
    ### RETRIEVE TAD PVALUES
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(file.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    tad_pval <- sort(tad_pval, decreasing=FALSE)
    
    if(nTopTADs == "all") {
      topTADs <- names(tad_pval)
    } else {
      topTADs <- names(tad_pval[1:nTopTADs])
    }
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    stopifnot(length(pipelineGenes) > 0)
    stopifnot(pipelineGenes %in% gene2tad_DT$entrezID)
    
    
    curr_g2t_DT <- gene2tad_DT[gene2tad_DT$entrezID %in% pipelineGenes &
                                 gene2tad_DT$region %in% names(tad_pval),
                                 ]
    stopifnot(nrow(curr_g2t_DT) > 0)
    nGenesByTAD <- setNames(as.numeric(table(curr_g2t_DT$region)), as.character(names(table(curr_g2t_DT$region))))
    
    stopifnot(topTADs %in% curr_g2t_DT$region)
    
    topTADs_nbrGenes <- nGenesByTAD[names(nGenesByTAD) %in% topTADs]
    stopifnot(length(topTADs_nbrGenes) == length(topTADs))  
    stopifnot(topTADs %in% names(topTADs_nbrGenes))  
      
    # random sampling of TAD size
    xx <- foreach(i = 1:nRandom) %dopar% {
      randSamp <- sample(nGenesByTAD, nTopTADs)
      names(randSamp) <- paste0("rand", i, "_TAD", 1:nTopTADs)
      randSamp
    }
      
    
    obsDT <- data.frame(
      dataset = curr_ds, 
      region = topTADs,
      nbrGenes = topTADs_nbrGenes[topTADs],
      tadType = "observed",
      stringsAsFactors=FALSE
    )
    
    randomDT <- data.frame(
      dataset = curr_ds,
      region = unlist(lapply(xx, function(x) names(x))),
      nbrGenes = unlist(lapply(xx, function(x) x)),
      tadType = "random",
      stringsAsFactors = FALSE
    )

    rbind(obsDT, randomDT)
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

# stop("-- ok -- \n")


############################################################################################################
############################################################################################################

all_ds_aucFCC <- foreach(curr_ds = unique(all_ds_DT$dataset), .combine='c') %dopar% {
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
names(all_ds_aucFCC) <- unique(all_ds_DT$dataset)
all_ds_aucFCC <- sort(all_ds_aucFCC, decreasing = TRUE)

topTADs <- names(all_ds_aucFCC[1:nToPlot])
botTADs <- names(rev(all_ds_aucFCC)[1:nToPlot])

############################################################################################################
############################################################################################################ all datasets
############################################################################################################

plotDT_m <- all_ds_DT

plotDT_m$dataset <- factor(plotDT_m$dataset, levels = names(all_ds_aucFCC))


plotTit <- paste0("#genes/TAD - obs./random")
mySub <- paste0("# top TADs = ", nTopTADs, "; # random = ", nRandom)

p_var <- ggplot(plotDT_m, aes(x = dataset, y = nbrGenes, fill = tadType)) +
  ggtitle(plotTit, subtitle = mySub)+
  geom_boxplot()+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("# genes/TAD"),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2"), unique(plotDT_m$tadType)),
                    labels = setNames(unique(plotDT_m$tadType), unique(plotDT_m$tadType)))+
  labs(fill = "")+
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size=10),
    panel.grid = element_blank(),
    axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=8),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    # axis.ticks.x = element_blank(),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank()
  )
if(SSHFS) p_var
outFile <- file.path(outFold, paste0("nbr_genes_by_TAD", "_boxplot_allDS.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


my_comparisons <- list(unique(plotDT_m$tadType))
p_box <- ggviolin(plotDT_m, x = "tadType", y = "nbrGenes", 
                  fill = "tadType",
                  palette = c("#00AFBB", "#FC4E07"),
                  title = plotTit,
                  xlab="",
                  ylab = paste0("# genes/TAD"),
                  # sub=mySub,
                  legend.title = "",
                  add = "boxplot", add.params = list(fill = "white"))+
  # stat_compare_means(comparisons = my_comparisons, 
  #                    # position = "bottomleft",
  #                    label = "p.signif")+ # Add significance levels
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, 
                     method="wilcox",
                     # position = "bottomleft",
                     label = "p.format")+ # Add significance levels
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

outFile <- file.path(outFold, paste0("nbr_genes_by_TAD", "_boxplot_byType_allDS.", plotType))
ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



############################################################################################################
############################################################################################################ top and bottom only
############################################################################################################


plotDT_m <- all_ds_DT[all_ds_DT$dataset %in% c(topTADs, botTADs),]

plotDT_m$dataset <- factor(plotDT_m$dataset, levels = c(topTADs, botTADs))

plotTit <- paste0("#genes/TAD - obs./random (", nToPlot, " top and bot. datasets)")
mySub <- paste0("# top TADs = ", nTopTADs, "; # random = ", nRandom)

p_var <- ggplot(plotDT_m, aes(x = dataset, y = nbrGenes, fill = tadType)) +
  ggtitle(plotTit, subtitle = mySub)+
  geom_boxplot()+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("# genes/TAD"),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2"), unique(plotDT_m$tadType)),
                    labels = setNames(unique(plotDT_m$tadType), unique(plotDT_m$tadType)))+
  labs(fill = "")+
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 2, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size=10),
    panel.grid = element_blank(),
    axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=8),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    # axis.ticks.x = element_blank(),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank()
  )
if(SSHFS) p_var
outFile <- file.path(outFold, paste0("nbr_genes_by_TAD", "_boxplot_topBotDS.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


plotDT_m <- all_ds_DT[all_ds_DT$dataset %in% c(topTADs),]
plotTit <- paste0("#genes/TAD - obs./random (", nToPlot, " top datasets)")

my_comparisons <- list(unique(plotDT_m$tadType))
p_box <- ggviolin(plotDT_m, x = "tadType", y = "nbrGenes", 
                  fill = "tadType",
                  palette = c("#00AFBB", "#FC4E07"),
                  title = plotTit,
                  xlab="",
                  ylab = paste0("# genes/TAD"),
                  # sub=mySub,
                  legend.title = "",
                  add = "boxplot", add.params = list(fill = "white"))+
  # stat_compare_means(comparisons = my_comparisons, 
  #                    # position = "bottomleft",
  #                    label = "p.signif")+ # Add significance levels
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, 
                     method="wilcox",
                     # position = "bottomleft",
                     label = "p.format")+ # Add significance levels
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

outFile <- file.path(outFold, paste0("nbr_genes_by_TAD", "_boxplot_byType_topDS.", plotType))
ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


plotDT_m <- all_ds_DT[all_ds_DT$dataset %in% c(botTADs),]
plotTit <- paste0("#genes/TAD - obs./random (", nToPlot, " bot. datasets)")

my_comparisons <- list(unique(plotDT_m$tadType))
p_box <- ggviolin(plotDT_m, x = "tadType", y = "nbrGenes", 
                  fill = "tadType",
                  palette = c("#00AFBB", "#FC4E07"),
                  title = plotTit,
                  xlab="",
                  ylab = paste0("# genes/TAD"),
                  # sub=mySub,
                  legend.title = "",
                  add = "boxplot", add.params = list(fill = "white"))+
  # stat_compare_means(comparisons = my_comparisons, 
  #                    # position = "bottomleft",
  #                    label = "p.signif")+ # Add significance levels
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, 
                     method="wilcox",
                     # position = "bottomleft",
                     label = "p.format")+ # Add significance levels
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

outFile <- file.path(outFold, paste0("nbr_genes_by_TAD", "_boxplot_byType_botDS.", plotType))
ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


