startTime <- Sys.time()
cat(paste0("> Rscript test_davidEnrich.R\n"))



options(scipen=100)

suppressPackageStartupMessages(library(RDAVIDWebService, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(reshape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
# suppressPackageStartupMessages(library(AnnotationDbi, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
# suppressPackageStartupMessages(library(GO.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
# suppressPackageStartupMessages(library(topGO, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) #inducedGraph
# suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(VennDiagram, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

vdHeight <- 7
vdWidth <- 7

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")
source("analysis_utils.R")

nTopDS <- 5
rankVar <- "FCC"


args <- commandArgs(trailingOnly = T)
if(length(args) >= 1) {
  if(!is.na(as.numeric(args[1])))
    nTopDS <- as.numeric(args[1])
}
if(length(args) == 2) {
  rankVar <- args[2]
}
stopifnot(rankVar %in% c("FCC", "coexpr", "avg"))

nTopResultsToSave <- 10

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05
pvalSelectGO <- 0.05

padjVarGO <- "p.adjust" # p.adjust or qvalue ???
barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

# settings for enrichment analysis: clusterProfiler::enricher
# default values enricher function: (universe as background !) 
# -> UNIVERSE changed 21.11: for the DE genes, it is genes used in limma analysis, for TAD DE genes it is genes from the TADs used in TAD DA analysis
# pvalueCutoff = 0.05 
# pAdjustMethod = "BH"
# minGSSize = 10
# maxGSSize = 500
# qvalueCutoff = 0.2 
enricher_ontologyType <- "BP"
enricher_pvalueCutoff <- 1
enricher_pAdjustMethod <- "BH"
enricher_minGSSize <- 1
enricher_maxGSSize <- 500
enricher_qvalueCutoff <- 1

outFold <- file.path("DAVID", enricher_ontologyType,nTopDS,rankVar)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) >= nTopDS)

#all_ds <- all_ds[all_ds %in% run1_DS]


txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

GO_g <- makeGOGraph(ont = tolower(enricher_ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)

enricher_results_sortGOby <- "p.adjust"

# GO for BP nad MF [do not take c5_CC]
if(enricher_ontologyType == "BP" | enricher_ontologyType == "MF" | enricher_ontologyType == "BP_MF"){
  gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", paste0("c5.", tolower(enricher_ontologyType), ".v6.1.entrez.gmt"))
} else {
  stop(paste0(enricher_ontologyType, " is not a valid ontologyType\n"))
}
stopifnot(file.exists(gmtFile))


mySub <- paste0("(pval thresh: signif. genes = ", pvalSelect , "; signif. GO = ",pvalSelectGO, ")") 

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... gmtFile:\t", gmtFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher background genes:\t", "universe", "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pvalueCutoff:\t", enricher_pvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pAdjustMethod:\t", enricher_pAdjustMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher minGSSize:\t", enricher_minGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher maxGSSize:\t", enricher_maxGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher qvalueCutoff:\t", enricher_qvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher ontologyType:\t", enricher_ontologyType, "\n")
printAndLog(txt, logFile)
txt <- paste0("... pval thresh select genes:\t", pvalSelect, "\n")
printAndLog(txt, logFile)
txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:5]

c5_msigdb <- read.gmt(gmtFile)
stopifnot(!is.numeric(c5_msigdb$gene))

nGenesByGo <- setNames(as.numeric(table(c5_msigdb$ont)), 
                       as.character(names(table(c5_msigdb$ont))))


##########################################################################################
##########################################################################################

# cat("all_ds[60] = ", all_ds[60], "\n")

stopifnot(!is.na(all_ds))

all_ds ="TCGAcrc_msi_mss"
curr_ds=all_ds

all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  if(!file.exists(aucFCC_file)) cat("aucFCC_file = ", aucFCC_file, "\n")
  stopifnot(file.exists(aucFCC_file))
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  aucFCC
}
names(all_ds_aucFCC) <- all_ds

all_ds_aucCoexprDist <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
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
names(all_ds_aucCoexprDist) <- all_ds

#cat("HELLO1\n")

# topDS <- c("GSE84231_lhb_rhb",
#            "GSE86356_tibMD1_tibNorm"
#            )

maxNds <- min(nTopDS, length(all_ds_aucFCC))

if(rankVar == "FCC"){
  topDS <- names(sort(all_ds_aucFCC, decreasing=TRUE)[1:maxNds])  
} else if(rankVar == "coexpr"){
  topDS <- names(sort(all_ds_aucCoexprDist, decreasing=TRUE)[1:maxNds])
} else if(rankVar == "avg"){
  tmpFCC <- sort(rank(-all_ds_aucFCC, ties = "first" ))
  tmpCoexpr <- sort(rank(-all_ds_aucCoexprDist, ties = "first" ))
  stopifnot(names(tmpFCC) %in% names(tmpCoexpr))
  stopifnot(names(tmpCoexpr) %in% names(tmpFCC))
  tmpFCC <- tmpFCC[names(tmpFCC)]
  tmpCoexpr <- tmpCoexpr[names(tmpFCC)]
  stopifnot(!is.na(tmpFCC))
  stopifnot(!is.na(tmpCoexpr))
  meanScore <- 0.5*(tmpFCC + tmpCoexpr)
  stopifnot(names(meanScore) == names(tmpFCC))
  topDS <- names(sort(meanScore, decreasing=F)[1:maxNds])
} else {
  stop("-- error rankVar -- \n")
}

stopifnot(!is.na(topDS))

##########################################################################################
##########################################################################################

#all_ds <- all_ds[1:3]
curr_ds="TCGAcrc_msi_mss"


  all_ds_DT <- foreach(curr_ds = topDS, .combine='rbind') %dopar% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    printAndLog(txt, logFile)
    
    ### RETRIEVE TAD PVALUES
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(file.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    # -> change here pval selection
    selectTADs <- names(tad_pval[tad_pval <= pvalSelect])
    nSelectTADs <- length(selectTADs)
    
    nTADs <- length(tad_pval)
    
    ### RETRIEVE GENES USED FOR THE LIMMA ANALYSIS - FOR THE UNIVERSE
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    rnageneFile <- file.path(step0_fold, "rna_geneList.Rdata")
    stopifnot(file.exists(rnageneFile))
    rnaGenes <- eval(parse(text = load(rnageneFile)))
    
    ### RETRIEVE GENES PVALUES
    step1_fold <- file.path(dsFold, curr_ds, "1_runGeneDE")
    stopifnot(file.exists(step1_fold))
    toptableFile <- file.path(step1_fold, "DE_topTable.Rdata")
    stopifnot(file.exists(toptableFile))
    topTable_DT <- eval(parse(text = load(toptableFile)))
    stopifnot(!any(duplicated(topTable_DT$genes)))
    stopifnot(topTable_DT$genes == rownames(topTable_DT))
    gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    
    ### SELECT UNIVERSE FOR THE DE GENES
    # stopifnot(names(rnaGenes) %in% topTable_DT$genes)# not true CPM filter ?
    stopifnot(topTable_DT$genes %in% names(rnaGenes))
    # stopifnot(length(rnaGenes) == nrow(topTable_DT)) # not true CPM filter ?
    stopifnot(length(rnaGenes) >= nrow(topTable_DT))
    universe_rnaGenes <- as.character(rnaGenes)
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    stopifnot(names(pipelineGenes) %in% topTable_DT$genes)
    
    
    ### SELECT UNIVERSE FOR THE TAD DE GENES
    universe_tadGenes <- as.character(pipelineGenes)
    
    topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipelineGenes),]
    genes_entrez <- sapply(topTable_DT$genes, function(x) pipelineGenes[x])
    
    stopifnot(!is.na(genes_entrez))
    stopifnot(length(genes_entrez) == nrow(topTable_DT) )
    stopifnot(names(pipelineGenes) %in% names(gene_pval))
    # gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    entrez_pval <- setNames(topTable_DT$adj.P.Val, genes_entrez)
    
    nGenes <- length(entrez_pval)
    
    selectGenes <- names(entrez_pval[entrez_pval <= pvalSelect])
    stopifnot(selectGenes %in% pipelineGenes)
    nSelectGenes <- length(selectGenes)
    
    stopifnot(selectTADs %in% gene2tad_DT$region)
    
    if(nSelectTADs > 0){
      stopifnot(selectTADs %in% gene2tad_DT$region)
      selectTADs_genes <- gene2tad_DT$entrezID[gene2tad_DT$region %in% selectTADs]  
      selectTADs_genes <- selectTADs_genes[selectTADs_genes %in% pipelineGenes]
      stopifnot(selectTADs_genes %in% gene2tad_DT$entrezID)
      stopifnot(selectTADs_genes %in% pipelineGenes)
      nSelectTADs_genes <- length(selectTADs_genes)
      
    } else {
      nSelectTADs_genes <- 0
    }
    
    if(nSelectGenes > 0){
      stopifnot(selectGenes %in% gene2tad_DT$entrezID)
      stopifnot(selectGenes %in% pipelineGenes)
    }
    
    
    
    david<-DAVIDWebService$new(email="marie.zufferey.1@unil.ch")
    data(demoList1)
    result<-addList(david, selectGenes,
                       idType="ENTREZ_GENE_ID",
                        listName="demoList1", listType="Gene")
    result
    
    listType="Background"
    .  If required,
    the  user  can  select  which  annotation  category  to  use,  e.g.
    GOTERM_BP_ALL
    ,
    GOTERM_MF_ALL
    and
    GOTERM_CC_ALL
    .
    R> setAnnotationCategories(david, c("GOTERM_BP_ALL",
                                        + "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
    Now, that everything is in order, the user can get the different reports to use
    right away or to save into a file for future recall.  For example the Functional
    Annotation  Clustering  can  be  obtained  on-line  on
    termCluster
    object,  or  as
    termClusterReport1.tab
    file by invoking:
      R> termCluster<-getClusterReport(david, type="Term")
    R> getClusterReportFile(david, type="Term",
                            + fileName="termClusterReport1.tab"
                            
                            