startTime <- Sys.time()
cat(paste0("> Rscript coexpr_dist_v3_sortNoDup.R\n"))

options(scipen=100)

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggstatsplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

axisLabSize <- 12
legendSize <- 10
plotTitSize <- 14

mytheme <- theme(
  # top, right, bottom and left
  plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
  plot.title = element_text(hjust = 0.5, face = "bold", size=plotTitSize, vjust=1),
  plot.subtitle = element_text(hjust = 0.5, face = "bold", size=plotTitSize-2, vjust=1),
  panel.background = element_rect(fill = "white", colour = NA), 
  panel.border = element_rect(fill = NA, colour = "grey20"), 
  panel.grid.major = element_line(colour = "grey92"), 
  panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
  strip.background = element_rect(fill = "grey85", colour = "grey20"), 
  #legend.key = element_rect(fill = "white", colour = NA), 
  axis.line.x = element_line(size = .3, color = "black"),
  axis.line.y = element_line(size = .3, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=axisLabSize),
  axis.text.x = element_text(color="black", hjust=0.5,vjust = 1, size=axisLabSize),
  axis.title.y = element_text(color="black", size=axisLabSize+1),
  axis.title.x = element_text(color="black", size=axisLabSize+1),
  legend.text = element_text(size=legendSize),
  legend.key.height = unit(1.5,"cm"),
  legend.key = element_blank()
)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

### HARD CODED
caller <- "TopDom"
corMethod <- "pearson"
buildTable <- TRUE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"

scatterFontSizeLabel <- 14
scatterFontSizeTitle <- 12

# UPDATE 30.06.2018:
# -> check that always $gene1 < $gene2 before left_join !!!

### RETRIEVE FROM COMMAND LINE
# Rscript coexpr_dist_v3.R

# top-ranking:
# Rscript coexpr_dist_v3_sortNoDup.R TCGAcrc_msi_mss hgnc # running

# Rscript coexpr_dist_v3_sortNoDup.R <dataset>
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 1)
  curr_dataset <- args[1]
} else{
  txt <- paste0("> Default parameters:\n")
  curr_dataset <- "TCGAcrc_msi_mss"
}

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "COEXPR_DIST_SAME_STRAND",  paste0(curr_dataset))
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("coexpr_dist_v3_logFile.txt"))  
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

mycols <- c("same TAD" ="darkorange1" , "diff. TAD"="darkslateblue",  "same Strand + same TAD"="violetred1", "same Strand + diff. TAD" = "lightskyblue")

plotType <- "svg"
# myHeight <- ifelse(plotType == "png", 400, 7)
# myWidth <- ifelse(plotType == "png", 600, 10)
myHeight <- ifelse(plotType == "png", 200, 5)
myWidth <- ifelse(plotType == "png", 350, 6)

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
  # UPDATE 30.06.2018
  stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
  
  cat(paste0("... load TAD data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_TAD_SORTNODUP/all_TAD_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_TAD_pairs)
  nrow(all_TAD_pairs)
  all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
  all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
  # UPDATE 30.06.2018
  stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
  
  cat(paste0("... load COEXPR data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_COEXPR_SORTNODUP"),  paste0(curr_dataset, "_", corMethod), "coexprDT.Rdata"))
  cat(paste0(Sys.time(), "\n"))
  head(coexprDT)
  nrow(coexprDT)
  coexprDT$gene1 <- as.character(coexprDT$gene1)
  coexprDT$gene2 <- as.character(coexprDT$gene2)
  all_TAD_pairs$gene2
  # UPDATE 30.06.2018
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
  head(dataset_TAD_pairs)
  nrow(dataset_TAD_pairs)
  
  
  # START MERGING DATA 
  
  cat(paste0("... merge DIST - TAD data\t", Sys.time(), "\t"))
  dataset_dist_TAD_DT <- left_join(dataset_dist_pairs_limit, dataset_TAD_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  dataset_dist_TAD_DT$sameTAD <- ifelse(is.na(dataset_dist_TAD_DT$region), 0, 1)
  
}



if(buildTable){
  
  
  cat(paste0("... load STRAND data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_SAME_STRAND_SORTNODUP", "all_strand_pairs.Rdata"))
  cat(paste0(Sys.time(), "\n"))
  head(all_strand_pairs)
  nrow(all_strand_pairs)
  all_strand_pairs$gene1 <- as.character(all_strand_pairs$gene1)
  all_strand_pairs$gene2 <- as.character(all_strand_pairs$gene2)
  stopifnot(all_strand_pairs$gene1 < all_strand_pairs$gene2)
  
  dataset_strand_pairs <- all_strand_pairs[all_strand_pairs$gene1 %in% pipeline_geneList & 
                                             all_strand_pairs$gene2 %in% pipeline_geneList,]
  head(dataset_strand_pairs)
  
  cat(paste0("... merge STRAND data\t", Sys.time(), "\t"))
  
  dataset_dist_TAD_strand_DT <- left_join(dataset_dist_TAD_DT, dataset_strand_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  
  cat("dataset_dist_TAD_strand_DT:\n")
  head(dataset_dist_TAD_strand_DT)
  
  stopifnot(!is.na(dataset_dist_TAD_strand_DT$sameStrand))
  
  cat(paste0("... merge COEXPR data\t", Sys.time(), "\t"))
  
  dataset_dist_TAD_strand_coexpr_DT <- left_join(dataset_dist_TAD_strand_DT, coexprDT, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  cat("dataset_dist_TAD_strand_coexpr_DT:\n")
  head(dataset_dist_TAD_strand_coexpr_DT)
  
  allData_dt <- dataset_dist_TAD_strand_coexpr_DT
  allData_dt$region <- NULL
  allData_dt <- na.omit(allData_dt)
  
  outFile <-file.path(outFold, paste0("allData_dt.Rdata"))
  save(allData_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))    
  
} else{
  outFile <-file.path(outFold, paste0("allData_dt.Rdata"))
  load(outFile)
  
}

nrow(allData_dt)
allData_dt$dist_kb <- allData_dt$dist/1000

allData_dt$curve1 <-  ifelse(allData_dt$sameTAD == "0", "diff. TAD", "same TAD")

allData_dt$curve2 <-  ifelse(allData_dt$sameStrand == "0", NA,
                             ifelse(allData_dt$sameTAD == "0", "same Strand + diff. TAD", "same Strand + same TAD"))

sameTAD_DT <- allData_dt[allData_dt$sameTAD == 1,c("gene1", "gene2", "coexpr", "dist_kb")]
sameTAD_DT <- na.omit(sameTAD_DT)
sameTAD_DT <- sameTAD_DT[order(sameTAD_DT$dist_kb),]
# sameTAD_DT$cumdist <- cumsum(sameTAD_DT$dist_kb)
sameTAD_DT$nPair <- 1:nrow(sameTAD_DT)
sameTAD_DT$label <- "same TAD"

diffTAD_DT <- allData_dt[allData_dt$sameTAD == 0,c("gene1", "gene2",  "coexpr","dist_kb")]
diffTAD_DT <- na.omit(diffTAD_DT)
diffTAD_DT <- diffTAD_DT[order(diffTAD_DT$dist_kb),]
# diffTAD_DT$cumdist <- cumsum(diffTAD_DT$dist_kb)
diffTAD_DT$nPair <- 1:nrow(diffTAD_DT)
diffTAD_DT$label <- "diff. TAD"

sameStrand_sameTAD_DT <- allData_dt[allData_dt$sameStrand == 1 & allData_dt$sameTAD == 1 ,c("gene1", "gene2", "coexpr", "dist_kb")]
sameStrand_sameTAD_DT <- na.omit(sameStrand_sameTAD_DT)
sameStrand_sameTAD_DT <- sameStrand_sameTAD_DT[order(sameStrand_sameTAD_DT$dist_kb),]
# sameStrand_sameTAD_DT$cumdist <- cumsum(sameStrand_sameTAD_DT$dist_kb)
sameStrand_sameTAD_DT$nPair <- 1:nrow(sameStrand_sameTAD_DT)
sameStrand_sameTAD_DT$label <- "same Strand + same TAD"

sameStrand_diffTAD_DT <- allData_dt[allData_dt$sameStrand == 1 & allData_dt$sameTAD == 0 ,c("gene1", "gene2",  "coexpr","dist_kb")]
sameStrand_diffTAD_DT <- na.omit(sameStrand_diffTAD_DT)
sameStrand_diffTAD_DT <- sameStrand_diffTAD_DT[order(sameStrand_diffTAD_DT$dist_kb),]
# sameStrand_diffTAD_DT$cumdist <- cumsum(sameStrand_diffTAD_DT$dist_kb)
sameStrand_diffTAD_DT$nPair <- 1:nrow(sameStrand_diffTAD_DT)
sameStrand_diffTAD_DT$label <- "same Strand + diff. TAD"

sum_4_curves_DT <- rbind(rbind(sameTAD_DT, diffTAD_DT),
                         rbind(sameStrand_sameTAD_DT,sameStrand_diffTAD_DT))

stopifnot(sum_4_curves_DT$dist_kb <= distLimit/1000)

sum_4_curves_DT$label <- factor(sum_4_curves_DT$label,
                                levels=names(mycols))

stopifnot(!any(is.na(sum_4_curves_DT$label)))

my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
my_xlab <- paste0("Distance between the 2 genes (kb)")
my_sub <- paste0(curr_dataset)



p_coexpr2curves <- my_ggscatterhist(
  sum_4_curves_DT[as.character(sum_4_curves_DT$label) %in% c("diff. TAD", "same TAD"),], 
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
outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameTADstrand_diffTADstrand_2curves_coexpr2_dist.", plotType))
ggsave(p_coexpr2curves, height = myHeight*1.2, width = myWidth*1.2, filename=outFile)
cat(paste0("... written: ", outFile, "\n"))


p_nbr <- ggplot(sum_4_curves_DT, aes(x=dist_kb, y = nPair, color=label))+
  ggtitle("Number of pairs along distance betw. genes") + 
  xlab(my_xlab)+
  ylab("Nbr of gene pairs")+
  geom_smooth(se=FALSE, method = fitMeth)+
  scale_color_manual(values=mycols, name="") + 
  mytheme

outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameTADstrand_diffTADstrand_4curves_nbr_dist.", plotType))
ggsave(p_nbr, height = myHeight, width = myWidth, filename=outFile)
cat(paste0("... written: ", outFile, "\n"))

p_nbr1 <- my_ggscatterhist(
  sum_4_curves_DT, 
  x = "dist_kb", 
  y = "nPair",
  xlab=my_xlab,
  ylab=my_ylab,
  title = paste0("Number of pairs along distance betw. genes"),
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
outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameTADstrand_diffTADstrand_4curves_nbr1_dist.", plotType))
ggsave(p_nbr1, height = myHeight*1.2, width = myWidth*1.2, filename=outFile)
cat(paste0("... written: ", outFile, "\n"))

p_coexpr1 <- ggplot(sum_4_curves_DT, aes(x=dist_kb, y = coexpr, color=label))+
  ggtitle("Coexpression along distance betw. genes", subtitle = my_sub) + 
  xlab(my_xlab)+
  ylab(my_ylab)+
  geom_smooth(se=FALSE, method = fitMeth)+
  scale_color_manual(values=mycols, name="") + 
  mytheme

outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameTADstrand_diffTADstrand_4curves_coexpr1_dist.", plotType))
ggsave(p_coexpr1, height = myHeight, width = myWidth, filename=outFile)
cat(paste0("... written: ", outFile, "\n"))


p_coexpr2 <- my_ggscatterhist(
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
  
  
  x_tit_font = c(scatterFontSizeTitle, "bold", "black"),
  x_text_font = c(scatterFontSizeLabel, "plain", "black"),
  y_tit_font = c(scatterFontSizeTitle, "plain", "black"),
  y_text_font = c(scatterFontSizeLabel, "plain", "black"),
  
  
  # global_margins_cm = c(0, 0, 0.25, 0.25)
  global_margins_cm = c(0.25, 0.25, 0.25, 0.25)
)
outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameTADstrand_diffTADstrand_4curves_coexpr2_dist.", plotType))
ggsave(p_coexpr2, height = myHeight*1.2, width = myWidth*1.2, filename=outFile)
cat(paste0("... written: ", outFile, "\n"))

p_coexpr <- ggplot() +
  ggtitle("Coexpression along distance betw. genes", subtitle = my_sub)+
  xlab(my_xlab)+
  ylab(my_ylab)+
  geom_smooth(data=allData_dt, aes(x=dist_kb, y=coexpr, color=curve1), 
              method=fitMeth, se=FALSE) +
  geom_smooth(data=allData_dt[allData_dt$sameStrand==1,], aes(x=dist_kb, y=coexpr, color=curve2), 
              method=fitMeth, se=FALSE)+
  scale_color_manual(values=mycols, name="") +
  mytheme

outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameTADstrand_diffTADstrand_4curves_coexpr_dist.", plotType))
ggsave(p_coexpr, height = myHeight, width = myWidth, filename=outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

