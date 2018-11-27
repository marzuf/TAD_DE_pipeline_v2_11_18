startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GO_top.R\n"))

# Rscript cmp_DE_TADs_genes_GO_top.R

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(reshape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
suppressPackageStartupMessages(library(AnnotationDbi, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
suppressPackageStartupMessages(library(GO.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
suppressPackageStartupMessages(library(topGO, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) #inducedGraph
suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(VennDiagram, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


# retrieved from ontologyIndex; needed for ontologyIndex::minimal_set()
data(go)

allGOs_AnnotationDBI <- keys(GO.db)
allGOs_ontologyIndex <- go$id

GOid_GOterm_DT_GOdb <- select(GO.db, keys=keys(GO.db), columns=c("GOID", "TERM"))
GOid_GOterm_DT_go <- data.frame(GOID = go$id, TERM=go$name, stringsAsFactors = F)
length(GOid_GOterm_DT_GOdb$TERM)
#44541
length(GOid_GOterm_DT_go$TERM)
# 45638
length(intersect(GOid_GOterm_DT_GOdb$TERM,GOid_GOterm_DT_go$TERM))
# 43442
intersectIDs <- intersect(GOid_GOterm_DT_GOdb$GOID,GOid_GOterm_DT_go$GOID)
length(intersectIDs)
# 43569
stopifnot(!any(duplicated(GOid_GOterm_DT_go$GOID)))
stopifnot(!any(duplicated(GOid_GOterm_DT_GOdb$GOID)))
# stopifnot(!any(duplicated(GOid_GOterm_DT_go$TERM)))
# stopifnot(!any(duplicated(GOid_GOterm_DT_GOdb$TERM)))

# they don't use the same terms (one uses abreviations in the terms, not the other)

GOid_GOterm_DT_cP <- clusterProfiler:::get_GO2TERM_table()
colnames(GOid_GOterm_DT_cP) <- c("GOID", "TERM")
stopifnot(!any(duplicated(GOid_GOterm_DT_cP$GOID)))
stopifnot(!any(duplicated(GOid_GOterm_DT_cP$TERM)))
GOterm_GOid_cP <- setNames(GOid_GOterm_DT_cP$GOID, GOid_GOterm_DT_cP$TERM)


length(allGOs_AnnotationDBI)
# 44541
length(allGOs_ontologyIndex)
# 45638
length(intersect(allGOs_ontologyIndex, allGOs_AnnotationDBI))
# 43569

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
args <- commandArgs(trailingOnly = T)
if(length(args) == 1) {
  if(!is.na(as.numeric(args[1])))
    nTopDS <- as.numeric(args[1])
}

nTopResultsToSave <- 10

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

pvalSelectGO <- 0.05

padjVarGO <- "p.adjust" # p.adjust or qvalue ???

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)


# settings for enrichment analysis: clusterProfiler::enricher
# default values enricher function: (universe as background !)
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

outFold <- file.path("CMP_TADs_GENES_GO_top", enricher_ontologyType, nTopDS)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) >= nTopDS)


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


mySub <- paste0("(select top percent = ", topPercent , "; signif. GO = ",pvalSelectGO, ")") 

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
txt <- paste0("... topPercent select genes:\t", topPercent, "\n")
printAndLog(txt, logFile)
txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
printAndLog(txt, logFile)

# all_ds <- all_ds[1:5]

c5_msigdb <- read.gmt(gmtFile)

##########################################################################################
##########################################################################################

all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
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

topDS <- names(sort(all_ds_aucFCC, decreasing=TRUE)[1:nTopDS])


##########################################################################################
##########################################################################################


#all_ds <- all_ds[1:3]
curr_ds="TCGAcrc_msi_mss"
if(buildTable){
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
    tad_pval_sort <- sort(tad_pval, decreasing=F)
    tad_pval_rank <- rank(tad_pval, ties.method="min")
    tad_pval_rank <- sort(tad_pval_rank)
    stopifnot( names(tad_pval_rank[tad_pval_rank==1])[1] == names(tad_pval_sort[1]) ) 
    
    nTopTADs <- round(topPercent * length(tad_pval_sort))
    stopifnot(nTopTADs > 0)
    
    nTADs <- length(tad_pval_sort)
    
    topTADs <- names(tad_pval_sort[1:nTopTADs])
    
    ### RETRIEVE GENES PVALUES
    step1_fold <- file.path(dsFold, curr_ds, "1_runGeneDE")
    stopifnot(file.exists(step1_fold))
    toptableFile <- file.path(step1_fold, "DE_topTable.Rdata")
    stopifnot(file.exists(toptableFile))
    topTable_DT <- eval(parse(text = load(toptableFile)))
    stopifnot(!any(duplicated(topTable_DT$genes)))
    stopifnot(topTable_DT$genes == rownames(topTable_DT))
    gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    stopifnot(names(pipelineGenes) %in% topTable_DT$genes)
    
    topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipelineGenes),]
    genes_entrez <- sapply(topTable_DT$genes, function(x) pipelineGenes[x])
    
    stopifnot(!is.na(genes_entrez))
    stopifnot(length(genes_entrez) == nrow(topTable_DT) )
    stopifnot(names(pipelineGenes) %in% names(gene_pval))
    # gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    entrez_pval <- setNames(topTable_DT$adj.P.Val, genes_entrez)
    entrez_pval_sort <- sort(entrez_pval, decreasing=F)
    entrez_pval_rank <- rank(entrez_pval, ties.method="min")
    
    stopifnot( names(entrez_pval_rank[entrez_pval_rank==1])[1] == names(entrez_pval_sort[1]) ) 
    stopifnot( names(entrez_pval_rank) %in% gene2tad_DT$entrezID )
    stopifnot( names(tad_pval_rank) %in% gene2tad_DT$region )
    
    nTopGenes <- round(topPercent * length(entrez_pval_sort))
    
    nGenes <- length(entrez_pval_sort)
    
    # take the top percent DE genes
    topGenes_manyAsPercent <- names(entrez_pval_sort[1:nTopGenes])
    nTopGenes_manyAsPercent <- length(topGenes_manyAsPercent)
    
    stopifnot(nTopGenes > 0)
    stopifnot(topTADs %in% gene2tad_DT$region)
    
    topTADs_genes <- gene2tad_DT$entrezID[gene2tad_DT$region %in% topTADs]  
    topTADs_genes <- topTADs_genes[topTADs_genes %in% pipelineGenes]
    nTopTADs_genes <- length(topTADs_genes)
    stopifnot(nTopTADs_genes > 0)
    stopifnot(nTopTADs_genes <= length(entrez_pval_sort) )
    
    # take as many top genes as genes from topTADs  
    topGenes_manyAsTopTADs <- names(entrez_pval_sort[1:length(topTADs_genes)])
    nTopGenes_manyAsTopTADs <- length(topGenes_manyAsTopTADs)
    
    stopifnot(length(topGenes_manyAsTopTADs) == length(topTADs_genes) )
    stopifnot(topTADs_genes %in% gene2tad_DT$entrezID)
    stopifnot(topTADs_genes %in% pipelineGenes)
    stopifnot(topGenes_manyAsPercent %in% gene2tad_DT$entrezID)
    stopifnot(topGenes_manyAsPercent %in% pipelineGenes)
    stopifnot(topGenes_manyAsTopTADs %in% gene2tad_DT$entrezID)
    stopifnot(topGenes_manyAsTopTADs %in% pipelineGenes)
  
    
    ############################################################################################
    ############################################################################################ enricher - topTADs
    ############################################################################################
    
  # run GO analysis for these 3 sets of entrez genes: topTADs_genes, topGenes_manyAsPercent
  #***** 1) topTADs_genes
  topTADs_genes_enrich <- enricher(gene = topTADs_genes, 
                                      TERM2GENE=c5_msigdb,
                                      pvalueCutoff = enricher_pvalueCutoff, 
                                      pAdjustMethod = enricher_pAdjustMethod, 
                                      minGSSize = enricher_minGSSize, 
                                      maxGSSize = enricher_maxGSSize, 
                                      qvalueCutoff =enricher_qvalueCutoff)
  
  topTADs_genes_resultDT <- topTADs_genes_enrich@result
  topTADs_genes_resultDT <- topTADs_genes_resultDT[order(topTADs_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]

  if(nrow(topTADs_genes_resultDT) < nTopResultsToSave) {
    topTADs_genes_enrichSaveDT <- topTADs_genes_resultDT
  } else {
    topTADs_genes_enrichSaveDT <- topTADs_genes_resultDT[1:nTopResultsToSave,]
  }
  outFile <- file.path(outFold, paste0(curr_ds, "_topTADs_genes_enrichSaveDT.Rdata"))
  save(topTADs_genes_enrichSaveDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))

  
  # topTADs_genes_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topTADs_genes_resultDT)))
  # if(topTADs_genes_nTop > 0) {
  #   topTADs_genes_topGo  <- as.character(topTADs_genes_resultDT$ID[1:topTADs_genes_nTop])
  # } else {
  #   topTADs_genes_topGo  <- character(0)  
  # }
  # stopifnot(length(topTADs_genes_topGo) == topTADs_genes_nTop)
  
  topTADs_genes_signifGOterm <- as.character(topTADs_genes_resultDT$ID[topTADs_genes_resultDT[, paste0(padjVarGO)] <= pvalSelectGO])
  nTopTADs_genes_signifGOterm <- length(topTADs_genes_signifGOterm)
  
  topTADs_genes_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", topTADs_genes_signifGOterm)))
  
  length(topTADs_genes_signifGOterm_tmp)
  # 58
  sum(topTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
  # 46
  sum(topTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
  # 46
  sum(topTADs_genes_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
  # 46
  
  topTADs_genes_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% topTADs_genes_signifGOterm_tmp]
  nTopTADs_genes_signifGOid <- length(topTADs_genes_signifGOid)
  
  # if(nTopTADs_genes_signifGOterm > 0)
  #   stopifnot(nTopTADs_genes_signifGOid > 0) # not always TRUE (task 3)
  

  ############################################################################################
  ############################################################################################ enricher - manyAsPercent
  ############################################################################################
  
#***** 2) topGenes_manyAsPercent
  topGenes_manyAsPercent_enrich <- enricher(gene = topGenes_manyAsPercent, 
                                 TERM2GENE=c5_msigdb,
                                 pvalueCutoff = enricher_pvalueCutoff, 
                                 pAdjustMethod = enricher_pAdjustMethod, 
                                 minGSSize = enricher_minGSSize, 
                                 maxGSSize = enricher_maxGSSize, 
                                 qvalueCutoff =enricher_qvalueCutoff)
  
  topGenes_manyAsPercent_resultDT <- topGenes_manyAsPercent_enrich@result
  topGenes_manyAsPercent_resultDT <- topGenes_manyAsPercent_resultDT[order(topGenes_manyAsPercent_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]

  if(nrow(topGenes_manyAsPercent_resultDT) < nTopResultsToSave) {
    topGenes_manyAsPercent_enrichSaveDT <- topGenes_manyAsPercent_resultDT
  } else {
    topGenes_manyAsPercent_enrichSaveDT <- topGenes_manyAsPercent_resultDT[1:nTopResultsToSave,]
  }
  outFile <- file.path(outFold, paste0(curr_ds, "_topGenes_manyAsPercent_enrichSaveDT.Rdata"))
  save(topGenes_manyAsPercent_enrichSaveDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))

  
  # topGenes_manyAsPercent_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topGenes_manyAsPercent_resultDT)))
  # if(topGenes_manyAsPercent_nTop > 0) {
  #   topGenes_manyAsPercent_topGo  <- as.character(topGenes_manyAsPercent_resultDT$ID[1:topGenes_manyAsPercent_nTop])
  # } else {
  #   topGenes_manyAsPercent_topGo  <- character(0)  
  # }
  # stopifnot(length(topGenes_manyAsPercent_topGo) == topGenes_manyAsPercent_nTop)
  
  topGenes_manyAsPercent_signifGOterm <- as.character(topGenes_manyAsPercent_resultDT$ID[topGenes_manyAsPercent_resultDT[, paste0(padjVarGO)] <= pvalSelectGO])
  nTopGenes_manyAsPercent_signifGOterm <- length(topGenes_manyAsPercent_signifGOterm)
  
  topGenes_manyAsPercent_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", topGenes_manyAsPercent_signifGOterm)))
  
  length(topGenes_manyAsPercent_signifGOterm_tmp)
  # 5
  sum(topGenes_manyAsPercent_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
  # 5
  sum(topGenes_manyAsPercent_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
  # 5
  sum(topGenes_manyAsPercent_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
  # 5
  
  topGenes_manyAsPercent_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% topGenes_manyAsPercent_signifGOterm_tmp]
  nTopGenes_manyAsPercent_signifGOid <- length(topGenes_manyAsPercent_signifGOid)
  
  # if(nTopGenes_manyAsPercent_signifGOterm > 0)
  #   stopifnot(nTopGenes_manyAsPercent_signifGOid > 0) # not true task 22
  
  ############################################################################################
  ############################################################################################ enricher - manyAsTopTADs
  ############################################################################################
  
  #***** 3) topGenes_manyAsTopTADs
  topGenes_manyAsTopTADs_enrich <- enricher(gene = topGenes_manyAsTopTADs, 
                                            TERM2GENE=c5_msigdb,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
  
  topGenes_manyAsTopTADs_resultDT <- topGenes_manyAsTopTADs_enrich@result
  topGenes_manyAsTopTADs_resultDT <- topGenes_manyAsTopTADs_resultDT[order(topGenes_manyAsTopTADs_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]

  if(nrow(topGenes_manyAsTopTADs_resultDT) < nTopResultsToSave) {
    topGenes_manyAsTopTADs_enrichSaveDT <- topGenes_manyAsTopTADs_resultDT
  } else {
    topGenes_manyAsTopTADs_enrichSaveDT <- topGenes_manyAsTopTADs_resultDT[1:nTopResultsToSave,]
  }
  outFile <- file.path(outFold, paste0(curr_ds, "_topGenes_manyAsTopTADs_enrichSaveDT.Rdata"))
  save(topGenes_manyAsTopTADs_enrichSaveDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))

  
  # topGenes_manyAsTopTADs_nTop <- min(c(enricher_results_cmpNbrTop, nrow(topGenes_manyAsTopTADs_resultDT)))
  # if(topGenes_manyAsTopTADs_nTop > 0) {
  #   topGenes_manyAsTopTADs_topGo  <- as.character(topGenes_manyAsTopTADs_resultDT$ID[1:topGenes_manyAsTopTADs_nTop])
  # } else {
  #   topGenes_manyAsTopTADs_topGo  <- character(0)  
  # }
  # stopifnot(length(topGenes_manyAsTopTADs_topGo) == topGenes_manyAsTopTADs_nTop)
  
  topGenes_manyAsTopTADs_signifGOterm <- as.character(topGenes_manyAsTopTADs_resultDT$ID[topGenes_manyAsTopTADs_resultDT[,paste0(padjVarGO)] <= pvalSelectGO])
  nTopGenes_manyAsTopTADs_signifGOterm <- length(topGenes_manyAsTopTADs_signifGOterm)
  
  topGenes_manyAsTopTADs_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", topGenes_manyAsTopTADs_signifGOterm)))
  
  length(topGenes_manyAsTopTADs_signifGOterm_tmp)
  # 5
  sum(topGenes_manyAsTopTADs_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
  # 5
  sum(topGenes_manyAsTopTADs_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
  # 5
  sum(topGenes_manyAsTopTADs_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
  # 5
  
  topGenes_manyAsTopTADs_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% topGenes_manyAsTopTADs_signifGOterm_tmp]
  nTopGenes_manyAsTopTADs_signifGOid <- length(topGenes_manyAsTopTADs_signifGOid)
  
  # if(nTopGenes_manyAsPercent_signifGOterm > 0)
  #   stopifnot(nTopGenes_manyAsPercent_signifGOid > 0) # not true task 22
  
  ############################################################################################
  ############################################################################################ GRAPH - topTADs
  ############################################################################################

  ###### GO analysis
  
  if(!is.na(nTopTADs_genes_signifGOid) & nTopTADs_genes_signifGOid > 0) {
    
    topTADs_genes_signifGOmin <- minimal_set(go, topTADs_genes_signifGOid)
    nTopTADs_genes_signifGOmin <- length(topTADs_genes_signifGOmin)
    topTADs_genes_gNEL <- inducedGraph(dag = GO_g, startNodes = topTADs_genes_signifGOmin[topTADs_genes_signifGOmin %in% nodes(GO_g)])
    
    # !!! need to reverse the graph !!!
    topTADs_genes_gNEL_rev <- reverseArch(topTADs_genes_gNEL) # topGO
    topTADs_genes_gNEL_rev_noAll <- removeNode("all", topTADs_genes_gNEL_rev) # not needed ? -> do not want this root
    topTADs_genes_g <- igraph.from.graphNEL(graphNEL = topTADs_genes_gNEL_rev_noAll)
    
    stopifnot(is.connected(topTADs_genes_g))
    stopifnot(is.directed(topTADs_genes_g))
    stopifnot(is.dag(topTADs_genes_g))
    
    topTADs_genes_g_nNodes <- length(V(topTADs_genes_g))
    
    # Now we can find the root nodes. (either no neighbors or no incident edges)
    # topTADs_genes_g_rootIdx <- which(sapply(sapply(V(topTADs_genes_g), 
    #                                 function(x) neighbors(topTADs_genes_g, x, mode="in")), length) == 0)
    # stopifnot(length(topTADs_genes_g_rootIdx) == 1)
    # topTADs_genes_g_rootGO <- names(V(topTADs_genes_g)[topTADs_genes_g_rootIdx])
    # txt <- paste0("...... found root GO: ", topTADs_genes_g_rootGO, "\n")
    # printAndLog(txt, log_file)
    
    # density
    topTADs_genes_g_density <- edge_density(topTADs_genes_g)
    
    # diameter
    topTADs_genes_g_diameter <- diameter(topTADs_genes_g)
    
    # mean geodesic distance
    topTADs_genes_g_meanDist <- mean_distance(topTADs_genes_g, directed = FALSE) # directed TRUE or FALSE ???
    topTADs_genes_g_meanDistDir <- mean_distance(topTADs_genes_g, directed = TRUE) # directed TRUE or FALSE ???
  
    # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
    topTADs_genes_g_meanEccentricity <- mean(eccentricity(topTADs_genes_g))
    
    # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
    topTADs_genes_g_meanBetweenness <- mean(betweenness(topTADs_genes_g))
    
    
  } else {
    topTADs_genes_signifGOmin <- NA
    nTopTADs_genes_signifGOmin <- NA
    topTADs_genes_g_nNodes <- NA
    topTADs_genes_g_density <- NA
    topTADs_genes_g_diameter <- NA
    topTADs_genes_g_meanDist <- NA
    topTADs_genes_g_meanDistDir <- NA
    topTADs_genes_g_meanEccentricity <- NA
    topTADs_genes_g_meanBetweenness <- NA
  }
    
  ############################################################################################
  ############################################################################################ GRAPH - manyAsPercent
  ############################################################################################
      
  if(!is.na(nTopGenes_manyAsPercent_signifGOid) & nTopGenes_manyAsPercent_signifGOid > 0) {
    
    topGenes_manyAsPercent_signifGOmin <- minimal_set(go, topGenes_manyAsPercent_signifGOid)
    nTopGenes_manyAsPercent_signifGOmin <- length(topGenes_manyAsPercent_signifGOmin)
    topGenes_manyAsPercent_gNEL <- inducedGraph(dag = GO_g, startNodes = topGenes_manyAsPercent_signifGOmin[topGenes_manyAsPercent_signifGOmin %in% nodes(GO_g)])
    
    # !!! need to reverse the graph !!!
    topGenes_manyAsPercent_gNEL_rev <- reverseArch(topGenes_manyAsPercent_gNEL) # topGO
    topGenes_manyAsPercent_gNEL_rev_noAll <- removeNode("all", topGenes_manyAsPercent_gNEL_rev) # not needed ? -> do not want this root
    topGenes_manyAsPercent_g <- igraph.from.graphNEL(graphNEL = topGenes_manyAsPercent_gNEL_rev_noAll)
    
    stopifnot(is.connected(topGenes_manyAsPercent_g))
    stopifnot(is.directed(topGenes_manyAsPercent_g))
    stopifnot(is.dag(topGenes_manyAsPercent_g))
    
    topGenes_manyAsPercent_g_nNodes <- length(V(topGenes_manyAsPercent_g))
    
    # Now we can find the root nodes. (either no neighbors or no incident edges)
    # topGenes_manyAsPercent_g_rootIdx <- which(sapply(sapply(V(topGenes_manyAsPercent_g), 
    #                                 function(x) neighbors(topGenes_manyAsPercent_g, x, mode="in")), length) == 0)
    # stopifnot(length(topGenes_manyAsPercent_g_rootIdx) == 1)
    # topGenes_manyAsPercent_g_rootGO <- names(V(topGenes_manyAsPercent_g)[topGenes_manyAsPercent_g_rootIdx])
    # txt <- paste0("...... found root GO: ", topGenes_manyAsPercent_g_rootGO, "\n")
    # printAndLog(txt, log_file)
    
    # density
    topGenes_manyAsPercent_g_density <- edge_density(topGenes_manyAsPercent_g)
    
    # diameter
    topGenes_manyAsPercent_g_diameter <- diameter(topGenes_manyAsPercent_g)
    
    # mean geodesic distance
    topGenes_manyAsPercent_g_meanDist <- mean_distance(topGenes_manyAsPercent_g, directed = FALSE) # directed TRUE or FALSE ???
    topGenes_manyAsPercent_g_meanDistDir <- mean_distance(topGenes_manyAsPercent_g, directed = TRUE) # directed TRUE or FALSE ???
    
    # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
    topGenes_manyAsPercent_g_meanEccentricity <- mean(eccentricity(topGenes_manyAsPercent_g))
    
    # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
    topGenes_manyAsPercent_g_meanBetweenness <- mean(betweenness(topGenes_manyAsPercent_g))
    
      
  } else {
    topGenes_manyAsPercent_signifGOmin <- NA
    nTopGenes_manyAsPercent_signifGOmin <- NA
    topGenes_manyAsPercent_g_nNodes <- NA
    topGenes_manyAsPercent_g_density <- NA
    topGenes_manyAsPercent_g_diameter <- NA
    topGenes_manyAsPercent_g_meanDist <- NA
    topGenes_manyAsPercent_g_meanDistDir <- NA
    topGenes_manyAsPercent_g_meanEccentricity <- NA
    topGenes_manyAsPercent_g_meanBetweenness <- NA
  }
    
    
  ############################################################################################
  ############################################################################################ GRAPH - manyAsTopTADs
  ############################################################################################
  
  if(!is.na(nTopGenes_manyAsTopTADs_signifGOid) & nTopGenes_manyAsTopTADs_signifGOid > 0) {
    
    topGenes_manyAsTopTADs_signifGOmin <- minimal_set(go, topGenes_manyAsTopTADs_signifGOid)
    nTopGenes_manyAsTopTADs_signifGOmin <- length(topGenes_manyAsTopTADs_signifGOmin)
    topGenes_manyAsTopTADs_gNEL <- inducedGraph(dag = GO_g, startNodes = topGenes_manyAsTopTADs_signifGOmin[topGenes_manyAsTopTADs_signifGOmin %in% nodes(GO_g)])
    
    # !!! need to reverse the graph !!!
    topGenes_manyAsTopTADs_gNEL_rev <- reverseArch(topGenes_manyAsTopTADs_gNEL) # topGO
    topGenes_manyAsTopTADs_gNEL_rev_noAll <- removeNode("all", topGenes_manyAsTopTADs_gNEL_rev) # not needed ? -> do not want this root
    topGenes_manyAsTopTADs_g <- igraph.from.graphNEL(graphNEL = topGenes_manyAsTopTADs_gNEL_rev_noAll)
    
    stopifnot(is.connected(topGenes_manyAsTopTADs_g))
    stopifnot(is.directed(topGenes_manyAsTopTADs_g))
    stopifnot(is.dag(topGenes_manyAsTopTADs_g))
    
    topGenes_manyAsTopTADs_g_nNodes <- length(V(topGenes_manyAsTopTADs_g))
    
    # Now we can find the root nodes. (either no neighbors or no incident edges)
    # topGenes_manyAsTopTADs_g_rootIdx <- which(sapply(sapply(V(topGenes_manyAsTopTADs_g), 
    #                                 function(x) neighbors(topGenes_manyAsTopTADs_g, x, mode="in")), length) == 0)
    # stopifnot(length(topGenes_manyAsTopTADs_g_rootIdx) == 1)
    # topGenes_manyAsTopTADs_g_rootGO <- names(V(topGenes_manyAsTopTADs_g)[topGenes_manyAsTopTADs_g_rootIdx])
    # txt <- paste0("...... found root GO: ", topGenes_manyAsTopTADs_g_rootGO, "\n")
    # printAndLog(txt, log_file)
    
    # density
    topGenes_manyAsTopTADs_g_density <- edge_density(topGenes_manyAsTopTADs_g)
    
    # diameter
    topGenes_manyAsTopTADs_g_diameter <- diameter(topGenes_manyAsTopTADs_g)
    
    # mean geodesic distance
    topGenes_manyAsTopTADs_g_meanDist <- mean_distance(topGenes_manyAsTopTADs_g, directed = FALSE) # directed TRUE or FALSE ???
    topGenes_manyAsTopTADs_g_meanDistDir <- mean_distance(topGenes_manyAsTopTADs_g, directed = TRUE) # directed TRUE or FALSE ???
    
    # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
    topGenes_manyAsTopTADs_g_meanEccentricity <- mean(eccentricity(topGenes_manyAsTopTADs_g))
    
    # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
    topGenes_manyAsTopTADs_g_meanBetweenness <- mean(betweenness(topGenes_manyAsTopTADs_g))
    
    
  } else {
    topGenes_manyAsTopTADs_signifGOmin <- NA
    nTopGenes_manyAsTopTADs_signifGOmin <- NA
    topGenes_manyAsTopTADs_g_nNodes <- NA
    topGenes_manyAsTopTADs_g_density <- NA
    topGenes_manyAsTopTADs_g_diameter <- NA
    topGenes_manyAsTopTADs_g_meanDist <- NA
    topGenes_manyAsTopTADs_g_meanDistDir <- NA
    topGenes_manyAsTopTADs_g_meanEccentricity <- NA
    topGenes_manyAsTopTADs_g_meanBetweenness <- NA
  }
  
  
  nIntersectGenes_manyAsPercent <- length(intersect(na.omit(topGenes_manyAsPercent), na.omit(topTADs_genes)))
  nIntersectGenes_manyAsTopTADs <- length(intersect(na.omit(topGenes_manyAsTopTADs), na.omit(topTADs_genes)))
  
  nUnionGenes_manyAsPercent <- length(union(na.omit(topGenes_manyAsPercent), na.omit(topTADs_genes)))
  nUnionGenes_manyAsTopTADs <- length(union(na.omit(topGenes_manyAsTopTADs), na.omit(topTADs_genes)))
  
  nIntersectSignifGOid_manyAsPercent <- length(intersect(na.omit(topTADs_genes_signifGOid), na.omit(topGenes_manyAsPercent_signifGOid)))
  nIntersectSignifGOterm_manyAsPercent <- length(intersect(na.omit(topTADs_genes_signifGOterm), na.omit(topGenes_manyAsPercent_signifGOterm)))
  nIntersectSignifGOmin_manyAsPercent <- length(intersect(na.omit(topTADs_genes_signifGOmin), na.omit(topGenes_manyAsPercent_signifGOmin)))
  nIntersectSignifGOid_manyAsTopTADs <- length(intersect(na.omit(topTADs_genes_signifGOid), na.omit(topGenes_manyAsTopTADs_signifGOid)))
  nIntersectSignifGOterm_manyAsTopTADs <- length(intersect(na.omit(topTADs_genes_signifGOterm), na.omit(topGenes_manyAsTopTADs_signifGOterm)))
  nIntersectSignifGOmin_manyAsTopTADs <- length(intersect(na.omit(topTADs_genes_signifGOmin), na.omit(topGenes_manyAsTopTADs_signifGOmin)))
  
  intersectSignifGOidRatio_manyAsPercent <- length(intersect(na.omit(topTADs_genes_signifGOid), na.omit(topGenes_manyAsPercent_signifGOid)))/
    length(union(na.omit(topTADs_genes_signifGOid), na.omit(topGenes_manyAsPercent_signifGOid)))
  intersectSignifGOtermRatio_manyAsPercent <- length(intersect(na.omit(topTADs_genes_signifGOterm), na.omit(topGenes_manyAsPercent_signifGOterm)))/
    length(union(na.omit(topTADs_genes_signifGOterm), na.omit(topGenes_manyAsPercent_signifGOterm)))
  intersectSignifGOminRatio_manyAsPercent <- length(intersect(na.omit(topTADs_genes_signifGOmin), na.omit(topGenes_manyAsPercent_signifGOmin)))/
    length(union(na.omit(topTADs_genes_signifGOmin), na.omit(topGenes_manyAsPercent_signifGOmin)))
  intersectSignifGOidRatio_manyAsTopTADs <- length(intersect(na.omit(topTADs_genes_signifGOid), na.omit(topGenes_manyAsTopTADs_signifGOid)))/
    length(union(na.omit(topTADs_genes_signifGOid), na.omit(topGenes_manyAsTopTADs_signifGOid)))
  intersectSignifGOtermRatio_manyAsTopTADs <- length(intersect(na.omit(topTADs_genes_signifGOterm), na.omit(topGenes_manyAsTopTADs_signifGOterm)))/
    length(union(na.omit(topTADs_genes_signifGOterm), na.omit(topGenes_manyAsTopTADs_signifGOterm)))
  intersectSignifGOminRatio_manyAsTopTADs <- length(intersect(na.omit(topTADs_genes_signifGOmin), na.omit(topGenes_manyAsTopTADs_signifGOmin)))/
    length(union(na.omit(topTADs_genes_signifGOmin), na.omit(topGenes_manyAsTopTADs_signifGOmin)))
  
  stopifnot(na.omit(intersectSignifGOidRatio_manyAsPercent) >= 0 & na.omit(intersectSignifGOidRatio_manyAsPercent <= 1 ))
  stopifnot(na.omit(intersectSignifGOtermRatio_manyAsPercent) >= 0 & na.omit(intersectSignifGOtermRatio_manyAsPercent <=1 ))
  stopifnot(na.omit(intersectSignifGOminRatio_manyAsPercent) >= 0 & na.omit(intersectSignifGOminRatio_manyAsPercent <=1 ))
  stopifnot(na.omit(intersectSignifGOidRatio_manyAsTopTADs) >= 0 & na.omit(intersectSignifGOidRatio_manyAsTopTADs <= 1 ))
  stopifnot(na.omit(intersectSignifGOtermRatio_manyAsTopTADs) >= 0 & na.omit(intersectSignifGOtermRatio_manyAsTopTADs <=1 ))
  stopifnot(na.omit(intersectSignifGOminRatio_manyAsTopTADs) >= 0 & na.omit(intersectSignifGOminRatio_manyAsTopTADs <=1 ))
  
  topGenes_manyAsPercent_intersectRatio <- nIntersectGenes_manyAsPercent/nTopGenes_manyAsPercent
  topTADs_genes_intersectRatio_manyAsPercent <- nIntersectGenes_manyAsPercent/nTopTADs_genes
  topGenes_manyAsTopTADs_intersectRatio <- nIntersectGenes_manyAsTopTADs/nTopGenes_manyAsTopTADs
  topTADs_genes_intersectRatio_manyAsTopTADs <- nIntersectGenes_manyAsTopTADs/nTopTADs_genes
  
  intersectGenesRatio_manyAsPercent <- nIntersectGenes_manyAsPercent/nUnionGenes_manyAsPercent
  intersectGenesRatio_manyAsTopTADs <- nIntersectGenes_manyAsTopTADs/nUnionGenes_manyAsTopTADs
  
  stopifnot(na.omit(topGenes_manyAsPercent_intersectRatio) >= 0 & na.omit(topGenes_manyAsPercent_intersectRatio <=1 ))
  stopifnot(na.omit(topTADs_genes_intersectRatio_manyAsPercent) >= 0 & na.omit(topTADs_genes_intersectRatio_manyAsPercent <=1 ))
  stopifnot(na.omit(intersectGenesRatio_manyAsPercent) >= 0 & na.omit(intersectGenesRatio_manyAsPercent <=1 ))
  stopifnot(na.omit(topGenes_manyAsTopTADs_intersectRatio) >= 0 & na.omit(topGenes_manyAsTopTADs_intersectRatio <=1 ))
  stopifnot(na.omit(topTADs_genes_intersectRatio_manyAsTopTADs) >= 0 & na.omit(topTADs_genes_intersectRatio_manyAsTopTADs <=1 ))
  stopifnot(na.omit(intersectGenesRatio_manyAsTopTADs) >= 0 & na.omit(intersectGenesRatio_manyAsTopTADs <=1 ))
  
  nTopTADs_genes_signifGOidRatio <- nTopTADs_genes_signifGOid/nTopTADs_genes
  nTopTADs_genes_signifGOtermRatio <- nTopTADs_genes_signifGOterm/nTopTADs_genes
  nTopTADs_genes_signifGOminRatio <- nTopTADs_genes_signifGOmin/nTopTADs_genes
  
  nTopGenes_manyAsPercent_signifGOidRatio <- nTopGenes_manyAsPercent_signifGOid/nTopGenes_manyAsPercent
  nTopGenes_manyAsPercent_signifGOtermRatio <- nTopGenes_manyAsPercent_signifGOterm/nTopGenes_manyAsPercent
  nTopGenes_manyAsPercent_signifGOminRatio <- nTopGenes_manyAsPercent_signifGOmin/nTopGenes_manyAsPercent
  
  nTopGenes_manyAsTopTADs_signifGOidRatio <- nTopGenes_manyAsTopTADs_signifGOid/nTopGenes_manyAsTopTADs
  nTopGenes_manyAsTopTADs_signifGOtermRatio <- nTopGenes_manyAsTopTADs_signifGOterm/nTopGenes_manyAsTopTADs
  nTopGenes_manyAsTopTADs_signifGOminRatio <- nTopGenes_manyAsTopTADs_signifGOmin/nTopGenes_manyAsTopTADs
  
    
    data.frame(
      dataset = curr_ds, 
      
    nTADs = nTADs,
    nTopTADs = nTopTADs,
      
    nTopTADs_genes=nTopTADs_genes,
    nTopTADs_genes_signifGOid=nTopTADs_genes_signifGOid,
    nTopTADs_genes_signifGOterm=nTopTADs_genes_signifGOterm,
    nTopTADs_genes_signifGOmin=nTopTADs_genes_signifGOmin,
    
    
    topTADs_genes_g_density=topTADs_genes_g_density,
    topTADs_genes_g_diameter=topTADs_genes_g_diameter,
    topTADs_genes_g_meanDist=topTADs_genes_g_meanDist,
    topTADs_genes_g_meanDistDir=topTADs_genes_g_meanDistDir,
    topTADs_genes_g_meanEccentricity=topTADs_genes_g_meanEccentricity,
    topTADs_genes_g_meanBetweenness=topTADs_genes_g_meanBetweenness,
    
    nGenes = nGenes,
    
    nTopGenes_manyAsPercent=nTopGenes_manyAsPercent,
    nTopGenes_manyAsPercent_signifGOid=nTopGenes_manyAsPercent_signifGOid,
    nTopGenes_manyAsPercent_signifGOterm=nTopGenes_manyAsPercent_signifGOterm,
    nTopGenes_manyAsPercent_signifGOmin=nTopGenes_manyAsPercent_signifGOmin,
    
    topGenes_manyAsPercent_g_density=topGenes_manyAsPercent_g_density,
    topGenes_manyAsPercent_g_diameter=topGenes_manyAsPercent_g_diameter,
    topGenes_manyAsPercent_g_meanDist=topGenes_manyAsPercent_g_meanDist,
    topGenes_manyAsPercent_g_meanDistDir=topGenes_manyAsPercent_g_meanDistDir,
    topGenes_manyAsPercent_g_meanEccentricity=topGenes_manyAsPercent_g_meanEccentricity,
    topGenes_manyAsPercent_g_meanBetweenness=topGenes_manyAsPercent_g_meanBetweenness,
    
    
    nTopGenes_manyAsTopTADs=nTopGenes_manyAsTopTADs,
    nTopGenes_manyAsTopTADs_signifGOid=nTopGenes_manyAsTopTADs_signifGOid,
    nTopGenes_manyAsTopTADs_signifGOterm=nTopGenes_manyAsTopTADs_signifGOterm,
    nTopGenes_manyAsTopTADs_signifGOmin=nTopGenes_manyAsTopTADs_signifGOmin,
    
    topGenes_manyAsTopTADs_g_density=topGenes_manyAsTopTADs_g_density,
    topGenes_manyAsTopTADs_g_diameter=topGenes_manyAsTopTADs_g_diameter,
    topGenes_manyAsTopTADs_g_meanDist=topGenes_manyAsTopTADs_g_meanDist,
    topGenes_manyAsTopTADs_g_meanDistDir=topGenes_manyAsTopTADs_g_meanDistDir,
    topGenes_manyAsTopTADs_g_meanEccentricity =topGenes_manyAsTopTADs_g_meanEccentricity,
    topGenes_manyAsTopTADs_g_meanBetweenness=topGenes_manyAsTopTADs_g_meanBetweenness,
    
    nIntersectSignifGOid_manyAsPercent = nIntersectSignifGOid_manyAsPercent,
    nIntersectSignifGOterm_manyAsPercent = nIntersectSignifGOterm_manyAsPercent,
    nIntersectSignifGOmin_manyAsPercent = nIntersectSignifGOmin_manyAsPercent,
    
    nIntersectSignifGOid_manyAsTopTADs = nIntersectSignifGOid_manyAsTopTADs,
    nIntersectSignifGOterm_manyAsTopTADs = nIntersectSignifGOterm_manyAsTopTADs,
    nIntersectSignifGOmin_manyAsTopTADs = nIntersectSignifGOmin_manyAsTopTADs,
    
    nTopTADs_genes_signifGOidRatio = nTopTADs_genes_signifGOidRatio,
    nTopTADs_genes_signifGOtermRatio = nTopTADs_genes_signifGOtermRatio,
    nTopTADs_genes_signifGOminRatio = nTopTADs_genes_signifGOminRatio,
    
    nTopGenes_manyAsPercent_signifGOidRatio = nTopGenes_manyAsPercent_signifGOidRatio,
    nTopGenes_manyAsPercent_signifGOtermRatio = nTopGenes_manyAsPercent_signifGOtermRatio,
    nTopGenes_manyAsPercent_signifGOminRatio = nTopGenes_manyAsPercent_signifGOminRatio,
    
    nTopGenes_manyAsTopTADs_signifGOidRatio = nTopGenes_manyAsTopTADs_signifGOidRatio,
    nTopGenes_manyAsTopTADs_signifGOtermRatio = nTopGenes_manyAsTopTADs_signifGOtermRatio,
    nTopGenes_manyAsTopTADs_signifGOminRatio = nTopGenes_manyAsTopTADs_signifGOminRatio,
    
    nIntersectGenes_manyAsPercent = nIntersectGenes_manyAsPercent,
    nIntersectGenes_manyAsTopTADs = nIntersectGenes_manyAsTopTADs,
    
    nUnionGenes_manyAsPercent = nUnionGenes_manyAsPercent,
    nUnionGenes_manyAsTopTADs = nUnionGenes_manyAsTopTADs,
    
    intersectSignifGOidRatio_manyAsPercent = intersectSignifGOidRatio_manyAsPercent,
    intersectSignifGOtermRatio_manyAsPercent = intersectSignifGOtermRatio_manyAsPercent,
    intersectSignifGOminRatio_manyAsPercent = intersectSignifGOminRatio_manyAsPercent,
    intersectSignifGOidRatio_manyAsTopTADs = intersectSignifGOidRatio_manyAsTopTADs,
    intersectSignifGOtermRatio_manyAsTopTADs = intersectSignifGOtermRatio_manyAsTopTADs,
    intersectSignifGOminRatio_manyAsTopTADs = intersectSignifGOminRatio_manyAsTopTADs,
    
    topGenes_manyAsPercent_intersectRatio = topGenes_manyAsPercent_intersectRatio,
    topTADs_genes_intersectRatio_manyAsPercent = topTADs_genes_intersectRatio_manyAsPercent,
    topGenes_manyAsTopTADs_intersectRatio = topGenes_manyAsTopTADs_intersectRatio,
    topTADs_genes_intersectRatio_manyAsTopTADs = topTADs_genes_intersectRatio_manyAsTopTADs,
    
    intersectGenesRatio_manyAsPercent = intersectGenesRatio_manyAsPercent,
    intersectGenesRatio_manyAsTopTADs = intersectGenesRatio_manyAsTopTADs,
    
    
    stringsAsFactors=FALSE
    
    )
    
  } # end iterate over all datasets
  
  rownames(all_ds_DT) <- NULL
  all_ds_DT$dataset <- as.character(all_ds_DT$dataset)

  # all_ds_DT <- all_ds_DT[order(all_ds_DT$nIntersectSignifGO, all_ds_DT$nIntersectSignifGOmin),]
  
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  save(all_ds_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_ds_DT.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_DT <- eval(parse(text = load(outFile)))
}

# stop("---OK ---\n")


# [1] "dataset"                        "nTopTADs_genes"             
# [3] "nTopTADs_genes_signifGOid"   "nTopTADs_genes_signifGOterm"
# [5] "nTopTADs_genes_signifGOmin"  "topTADs_genes_g_density"    
# [9] "topTADs_genes_g_diameter"    "topTADs_genes_g_meanDist"   
# [11] "topTADs_genes_g_meanDistDir" "nTopGenes_manyAsPercent"                  
# [13] "nTopGenes_manyAsPercent_signifGOid"        "nTopGenes_manyAsPercent_signifGOterm"     
# [15] "nTopGenes_manyAsPercent_signifGOmin"       "topGenes_manyAsPercent_g_density"         
# [19] "topGenes_manyAsPercent_g_diameter"         "topGenes_manyAsPercent_g_meanDist"        
# [21] "topGenes_manyAsPercent_g_meanDistDir"      "nIntersectSignifGOid"          
# [23] "nIntersectSignifGOterm"         "nIntersectSignifGOmin"         

stopifnot(!duplicated(all_ds_DT$dataset))  

# all_ds_DT$nTopTADs_genes_signifGOidRatio <- all_ds_DT$nTopTADs_genes_signifGOid/all_ds_DT$nTopTADs_genes
# all_ds_DT$nTopTADs_genes_signifGOtermRatio <- all_ds_DT$nTopTADs_genes_signifGOterm/all_ds_DT$nTopTADs_genes
# all_ds_DT$nTopTADs_genes_signifGOminRatio <- all_ds_DT$nTopTADs_genes_signifGOmin/all_ds_DT$nTopTADs_genes
# 
# all_ds_DT$nTopGenes_manyAsPercent_signifGOidRatio <- all_ds_DT$nTopGenes_manyAsPercent_signifGOid/all_ds_DT$nTopGenes_manyAsPercent
# all_ds_DT$nTopGenes_manyAsPercent_signifGOtermRatio <- all_ds_DT$nTopGenes_manyAsPercent_signifGOterm/all_ds_DT$nTopGenes_manyAsPercent
# all_ds_DT$nTopGenes_manyAsPercent_signifGOminRatio <- all_ds_DT$nTopGenes_manyAsPercent_signifGOmin/all_ds_DT$nTopGenes_manyAsPercent
# 
# all_ds_DT$nTopGenes_manyAsTopTADs_signifGOidRatio <- all_ds_DT$nTopGenes_manyAsTopTADs_signifGOid/all_ds_DT$nTopGenes_manyAsTopTADs
# all_ds_DT$nTopGenes_manyAsTopTADs_signifGOtermRatio <- all_ds_DT$nTopGenes_manyAsTopTADs_signifGOterm/all_ds_DT$nTopGenes_manyAsTopTADs
# all_ds_DT$nTopGenes_manyAsTopTADs_signifGOminRatio <- all_ds_DT$nTopGenes_manyAsTopTADs_signifGOmin/all_ds_DT$nTopGenes_manyAsTopTADs


stopifnot(all_ds_DT$dataset %in% names(all_ds_aucFCC))
stopifnot(all_ds_DT$dataset %in% names(all_ds_aucCoexprDist))
all_ds_DT$aucFCC <- all_ds_aucFCC[all_ds_DT$dataset]
all_ds_DT$aucCoexprDist <- all_ds_aucCoexprDist[all_ds_DT$dataset]

outDT <- all_ds_DT
outDT <- outDT[order(outDT$aucFCC, decreasing=TRUE),]
write.table(outDT, file = file.path(outFold, "all_ds_DT.txt"), sep="\t", quote=F, col.names=T, row.names=F)

myGenSub <- paste0("topPercent = ", topPercent, "; pvalSelectGO = ", pvalSelectGO)


######################################################################################
###################################################################################### VENN DIAGRAM
######################################################################################
# grid.newpage()
# draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), # => this will draw 11-11-9
#                    lty = rep("blank", 2), 
#                    fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
#                    cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

curr_ds="TCGAcrc_msi_mss"
for(curr_ds in unique(all_ds_DT$dataset)) {
  geneType = "manyAsPercent"
  for(geneType in c("manyAsPercent", "manyAsTopTADs")){
    
    all_VDs <- list()
    
    ds_dt <- all_ds_DT[all_ds_DT$dataset == curr_ds,]
    stopifnot(nrow(ds_dt) == 1)
    
    all_VDs[["nGenes"]] <- draw_myVenn(var1=paste0("nTopGenes_", geneType),
                                       var2="nTopTADs_genes", 
                                       varIntersect = paste0("nIntersectGenes_", geneType),
                                       tit = curr_ds,
                                       # subTit = paste0(GO_aliases_common_top[gsub(paste0("Genes_", geneType), "", paste0("nTopGenes_", geneType))], " - ", geneType),
                                       subTit = paste0(GO_aliases_common_top[gsub(paste0("Genes_", geneType), "", paste0("nTopGenes_", geneType))]),
                                       dataDT = ds_dt)
    
    curr_var="signifGOid"
    for(curr_var in c("signifGOid", "signifGOterm", "signifGOmin")) {
      all_VDs[[paste0(curr_var)]] <- draw_myVenn(var1=paste0("nTopGenes_", geneType, "_", curr_var),
                                                 var2=paste0("nTopTADs_genes_", curr_var),
                                                 varIntersect = paste0("nIntersect", capitalize(curr_var), "_", geneType),
                                                 tit = curr_ds,
                                                 # subTit = paste0(GO_aliases_common_top[paste0(curr_var)], " - ", geneType),
                                                 subTit = paste0(GO_aliases_common_top[paste0(curr_var)]),
                                                 dataDT = ds_dt,
                                                 subPatt1 = paste0("_", curr_var),
                                                 subPatt2 = paste0("_", curr_var)
      )
    } # end iterating over variables
    titGrob <- textGrob(
                        paste0(curr_ds),
                        # paste0(curr_ds, " - ", geneType),
                        gp=gpar(fontsize=20, fontface="bold"))
    all_vds <- do.call(arrangeGrob, c(all_VDs, nrow=2))
    finalVD <- arrangeGrob(all_vds, ncol = 1, top = titGrob)
    
    outFile <- file.path(outFold, paste0("intersect_", geneType, "_", curr_ds, "_venn_diagram.", plotType))
    ggsave(finalVD, file = outFile, height = vdHeight, width = vdWidth)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  } # end iterating over manyAsPercent or manyAsTopTADs
  

  
} # end iterating over datasets

#cat("HELLO4\n")

# plot all datasets

######################################################################################
###################################################################################### BARPLOT INTERSECT
######################################################################################

for(geneType in c("manyAsTopTADs", "manyAsPercent")) {
  intersect_dt <- all_ds_DT[, c(which(colnames(all_ds_DT) == "dataset"), grep(paste0("intersect.+Ratio_", geneType), colnames(all_ds_DT)))]
  intersect_dt_m <- melt(intersect_dt, id="dataset")
  stopifnot(intersect_dt_m$variable %in% names(GO_aliases_top))
  intersect_dt_m$variable_leg <- GO_aliases_top[as.character(intersect_dt_m$variable)]
  intersect_dt_m$variable_leg <- gsub(paste0(" - intersect ", geneType), "", intersect_dt_m$variable_leg)
  
  plotTit_intersect <- paste0("Intersect ratio - ", geneType)
  mySub_intersect <- paste0("all datasets (n=", length(unique(all_ds_DT$dataset)), ")")
  
  p_intersect <- ggplot(intersect_dt_m, aes(x = dataset, y = value, fill = variable_leg)) +
    ggtitle(plotTit_intersect, subtitle = mySub_intersect)+
    geom_bar(stat="identity", position = "dodge")+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0("intersect ratio"),
                       breaks = scales::pretty_breaks(n = 10))+
    # scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2"), unique(intersect_dt_m$variable_leg)),
    #                   labels = setNames(unique(intersect_dt_m$variable_leg), unique(intersect_dt_m$variable_leg)))+
    labs(fill = "")+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      plot.margin = unit(c(1, 1, 2, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text(color="black", hjust=1,vjust = 0.5, angle = 90, size=10),
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
  if(SSHFS) p_intersect
  
  outFile <- file.path(outFold, paste0("all_ds_intersect_",geneType, "_barplot.", plotType))
  ggsave(plot = p_intersect, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
}

# stop("--ok--\n")

######################################################################################
###################################################################################### BARPLOT P VALUES - 1 TYPE SEPARATELY
######################################################################################


# CMP_TADs_GENES_GO_top/20/TCGAluad_luad_mutKRAS_topGenes_manyAsPercent_enrichSaveDT.Rdata
curr_ds="TCGAcrc_msi_mss"
curr_type="topTADs_genes"

for(curr_ds in unique(all_ds_DT$dataset)) {
  for(curr_type in c("topTADs_genes", "topGenes_manyAsPercent", "topGenes_manyAsTopTADs")) {
    
    resultDT_file <- file.path(outFold, paste0(curr_ds, "_",curr_type, "_enrichSaveDT.Rdata"))
    
    cat(paste0("resultDT_file = ", resultDT_file, "\n"))
    
    if(! file.exists(resultDT_file)) next
    
    stopifnot(file.exists(resultDT_file))
    
    
    resultDT <- eval(parse(text = load(resultDT_file)))    
    
    stopifnot(nrow(resultDT) == nTopResultsToSave)
    
    resultDT$log10_pval <- -log10(resultDT$qvalue)
    
    myTit <- paste0(curr_ds, " - ", curr_type)
    
    outFile <- file.path(outFold,paste0(curr_ds, "_", curr_type, "_top", nTopResultsToSave, "_GOpvals.", plotType))
    do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
    par(oma=c(10,1,1,1))
    barplot(resultDT$log10_pval, 
            main = myTit, 
            ylab = paste0("-log10 (", padjVarGO, ")"),
            names.arg = gsub("GO_", "", resultDT$Description), las=2,
            cex.names = 0.6
    )
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}


######################################################################################
###################################################################################### BARPLOT P VALUES - ALL TOGETHER - ITERATE OVER DATASETS - ADDED
######################################################################################

# CMP_TADs_GENES_GO_top/20/TCGAluad_luad_mutKRAS_topGenes_manyAsPercent_enrichSaveDT.Rdata
curr_ds="TCGAcrc_msi_mss"
curr_type="topTADs_genes"

barplotSubtit <- paste0("topPercent=", topPercent, "; pvalSelectGO=", pvalSelectGO)

all_ds_pval_DT <- foreach(curr_ds = unique(all_ds_DT$dataset), .combine='rbind') %dopar% {
  ds_pval_dt <- foreach(curr_type = c("topGenes_manyAsPercent", "topGenes_manyAsTopTADs", "topTADs_genes"), .combine='rbind') %do% {
    
    resultDT_file <- file.path(outFold, paste0(curr_ds, "_",curr_type, "_enrichSaveDT.Rdata"))
    cat(paste0("resultDT_file = ", resultDT_file, "\n"))
    
    if(!file.exists(resultDT_file)) {
      result_DT <- data.frame(
        ID = NULL,
        log10_pval = NULL,
        rank = NULL,
        geneType = NULL,
        stringsAsFactors=F
      )
      return(result_DT)
    }
    
    stopifnot(file.exists(resultDT_file))
    resultDT <- eval(parse(text = load(resultDT_file)))    
    stopifnot(nrow(resultDT) == nTopResultsToSave)
    # resultDT$log10_pval <- -log10(resultDT$qvalue)
    resultDT$log10_pval <- -log10(resultDT[,paste0(padjVarGO)])
    resultDT <- resultDT[order(resultDT$log10_pval, decreasing=T),]
    resultDT <- resultDT[, c("ID", "log10_pval")]
    stopifnot(nrow(resultDT) > 0)
    resultDT$rank <- 1:nrow(resultDT)
    resultDT$geneType <- curr_type
    rownames(resultDT) <- NULL
    resultDT
  }
  
  ds_pval_dt <- ds_pval_dt[order(ds_pval_dt$rank, ds_pval_dt$geneType),]
  ds_pval_dt$x_ID <- c(1:nrow(ds_pval_dt))
  ds_pval_dt$x_ID <- factor(ds_pval_dt$x_ID, levels = unique(as.character(ds_pval_dt$x_ID)))
  ds_pval_dt$x_labels <- gsub("_", " ", gsub("GO_", "", as.character(ds_pval_dt$ID)))
  
  p_DS <- ggplot(ds_pval_dt, aes(x = x_ID, y = log10_pval, fill = geneType)) +
    facet_grid(~rank, switch="x", scale="free")+
    geom_bar(stat="identity", position="dodge") +
    scale_x_discrete(name="", labels=ds_pval_dt$x_labels, breaks = ds_pval_dt$x_ID)   +
    scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
                       breaks = scales::pretty_breaks(n = 5))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(
      values = setNames(c("dodgerblue4", "darkorange2", "chartreuse4"), unique(ds_pval_dt$geneType)),
      labels = setNames(unique(ds_pval_dt$geneType), unique(ds_pval_dt$geneType)) )+
    labs(fill  = "") +
    ggtitle(label = paste0("Top GO pval - ", curr_ds), subtitle = barplotSubtit) +
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      strip.text.x = element_text(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      #    axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    ) #+
  # geom_hline(yintercept = 1, linetype = 2)
  
  if(SSHFS) p_DS
  
  outFile <- file.path(outFold,paste0(curr_ds, "_", curr_type, "_top", nTopResultsToSave, "_GOpvals_topGenes_topTADs.", plotType))
  ggsave(p_DS, filename = outFile, height = myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  ds_pval_dt$dataset <- curr_ds
  ds_pval_dt
  
} # end iterating over datasets

rownames(all_ds_pval_DT) <- NULL
outFile <- file.path(outFold, "all_ds_pval_DT.Rdata")
save(all_ds_pval_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
###################################################################################### BARPLOT COMPARISON TOP GO PVALUES - ALL DATASETS - ADDED
######################################################################################

all_ds_pval_DT <- eval(parse(text = load(outFile)))

all_ds_pval_DT$rank <- factor(all_ds_pval_DT$rank, levels = as.character(sort(unique(as.numeric(as.character(all_ds_pval_DT$rank))))))

p_all <- ggplot(all_ds_pval_DT, aes(x = rank, y = log10_pval, fill = geneType, col = geneType)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter() +
  scale_x_discrete(name="GO rank")+
  scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
                     breaks = scales::pretty_breaks(n = 5))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
  # coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2", "chartreuse4"), unique(all_ds_pval_DT$geneType)),
                    labels = setNames(unique(all_ds_pval_DT$geneType), unique(all_ds_pval_DT$geneType)))+
  scale_colour_manual(values = setNames(c("dodgerblue4", "darkorange2", "chartreuse4"), unique(all_ds_pval_DT$geneType)),
                      labels = setNames(unique(all_ds_pval_DT$geneType), unique(all_ds_pval_DT$geneType)), guide = F)+
  labs(fill  = "") +
  ggtitle(label = paste0("Top GO pval - all DS (n=", length(unique(all_ds_pval_DT$dataset)), ")"), subtitle = barplotSubtit) +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size=10),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(),
    axis.text.x = element_text( hjust=0.5,vjust = 1, size=12),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    #    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank()
    # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
  ) #+
# geom_hline(yintercept = 1, linetype = 2)

if(SSHFS) p_all

outFile <- file.path(outFold,paste0("allDS_typeComp", "_top", nTopResultsToSave, "_GOpvals_boxplot.", plotType))
ggsave(p_all, filename = outFile, height = myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

p_all <- ggplot(all_ds_pval_DT, aes(x = rank, y = log10_pval, fill = geneType)) +
  geom_boxplot(outlier.shape=NA) +
  #geom_jitter() +
  scale_x_discrete(name="GO rank")+
  scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
                     breaks = scales::pretty_breaks(n = 5))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
  # coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2", "chartreuse4"), unique(all_ds_pval_DT$geneType)),
                    labels = setNames(unique(all_ds_pval_DT$geneType), unique(all_ds_pval_DT$geneType)))+
  labs(fill  = "") +
  ggtitle(label = paste0("Top GO pval - all DS (n=", length(unique(all_ds_pval_DT$dataset)), ")"), subtitle = barplotSubtit) +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size=10),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(),
    axis.text.x = element_text( hjust=0.5,vjust = 1, size=12),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    #    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank()
    # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
  ) #+
# geom_hline(yintercept = 1, linetype = 2)

if(SSHFS) p_all

outFile <- file.path(outFold,paste0("allDS_typeComp", "_top", nTopResultsToSave, "_GOpvals_boxplot_nojitter.", plotType))
ggsave(p_all, filename = outFile, height = myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
###################################################################################### SCATTERPLOT AUC FCC VS. AUC COEXPR DIST.
######################################################################################
# scatterplot vs. aucFCC and aucCoexprDist


all_vars <- colnames(all_ds_DT)
all_vars <- all_vars[all_vars != "dataset"]

ref_vars <- c("aucFCC", "aucCoexprDist")

stopifnot(all_vars %in% names(GO_offSets_top))
stopifnot(ref_vars %in% names(GO_offSets_top))
stopifnot(all_vars %in% names(GO_legPos_top))
stopifnot(ref_vars %in% names(GO_legPos_top))


ref_var=ref_vars[1]
curr_var=all_vars[1]

for(ref_var in ref_vars){
  for(curr_var in all_vars) {
    stopifnot(ref_var %in% colnames(all_ds_DT))
    stopifnot(curr_var %in% colnames(all_ds_DT))
  }
}

for(ref_var in ref_vars){
  for(curr_var in all_vars) {
    
    if(ref_var == curr_var) next
    
    if(grepl("topTADs_genes", curr_var) | grepl("TopTADs_genes", curr_var)) {
      mysub <- paste0("topTADs_genes - ", myGenSub)
    } else if(grepl("topGenes_manyAsPercent", curr_var) | grepl("TopGenes_manyAsPercent", curr_var) | grepl("_manyAsPercent", curr_var)) {
      mysub <- paste0("topGenes_manyAsPercent - ", myGenSub)
    } else if(grepl("topGenes_manyAsTopTADs", curr_var) | grepl("TopGenes_manyAsTopTADs", curr_var) | grepl("_manyAsTopTADs", curr_var)) {
      mysub <- paste0("topGenes_manyAsTopTADs - ", myGenSub)
    } else {
      mysub <- paste0(myGenSub)
    }
    
    # myylab <- unique(gsub("topTADs_genes_", "", 
    #                         gsub("topGenes_manyAsPercent_", "", 
    #                              gsub("topGenes_manyAsTopTADs_", "", 
    #                              gsub("nTopGenes_manyAsPercent_", "",
    #                                   gsub("nTopGenes_manyAsTopTADs_", "",
    #                                   gsub("nTopTADs_genes_", "",  curr_var)))))))
    # 
    myylab <- GO_aliases_top[curr_var]
    myxlab <- GO_aliases_top[ref_var]
    
    
    
    mymain <- paste0(myylab, " vs. ", ref_var)
    
    myx <- all_ds_DT[,ref_var]
    myy <- all_ds_DT[,curr_var]
    
    if(all(is.na(myx)) | all(is.na(myy))) next
    curr_colors <- dataset_proc_colors[as.character(all_ds_DT$dataset)]
    
    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(x=myx,
         y=myy,
         main = mymain,
         xlab =  myxlab,
         ylab = myylab,
         pch=16,cex=0.7)
    add_curv_fit(x=myx,
                 y=myy, withR2 = F, lty=2)
    addCorr(x=myx,
            y=myy, legPos = "topleft", bty="n")
    text(x = myx, 
         y = myy, 
         labels = all_ds_DT[,"dataset"], 
         col=curr_colors,
         cex=0.7)
    mtext(mysub,side=3)
    
    
    legend(
      GO_legPos_top[curr_var], 
      #"bottomright",
      legend=names(my_colors),
      lty=-1,
      pch=16,
      col = my_colors,
      lwd = 5,
      bty="n",
      cex = 0.6)
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

########################################## ########################################## 

all_cmp_var <- c("topGenes_manyAsPercent",  "topGenes_manyAsTopTADs", "topTADs_genes")

all_vars_toplot <- unique(gsub("topTADs_genes_", "", 
                      gsub("topGenes_manyAsPercent_", "", 
                           gsub("topGenes_manyAsTopTADs_", "", 
                           gsub("nTopGenes_manyAsPercent_", "",
                                gsub("nTopGenes_manyAsTopTADs_", "",
                                gsub("nTopTADs_genes_", "",  colnames(all_ds_DT))))))))


all_vars_toplot <- all_vars_toplot[! all_vars_toplot %in% c("dataset", "aucFCC", "aucCoexprDist",
                                                            
                                                            "nTopGenes_manyAsPercent",  "nTopGenes_manyAsTopTADs", 
                                                            "nTopTADs_genes",
                                                            
                                                            "nIntersectSignifGOterm_manyAsPercent","nIntersectSignifGOterm_manyAsTopTADs",
                                                            "nIntersectSignifGOid_manyAsPercent", "nIntersectSignifGOid_manyAsTopTADs", 
                                                            "nIntersectSignifGOmin_manyAsPercent","nIntersectSignifGOmin_manyAsTopTADs",
                                                            
                                                            "intersectSignifGOidRatio_manyAsPercent", "intersectSignifGOidRatio_manyAsTopTADs",
                                                            "intersectSignifGOminRatio_manyAsPercent", "intersectSignifGOminRatio_manyAsTopTADs",
                                                            "intersectSignifGOtermRatio_manyAsPercent", "intersectSignifGOtermRatio_manyAsTopTADs",
                                                            
                                                            "intersectGenesRatio_manyAsPercent", "intersectGenesRatio_manyAsTopTADs",
                                                            
                                                            "intersectRatio_manyAsPercent", "intersectRatio_manyAsTopTADs",
                                                            
                                                            "intersectRatio",
                                                            
                                                            "topGenes_manyAsPercent_intersectRatio", "topGenes_manyAsTopTADs_intersectRatio",
                                                            
                                                            "nIntersectGenes_manyAsPercent", "nIntersectGenes_manyAsTopTADs",
                                                            
                                                            "nUnionGenes_manyAsPercent","nUnionGenes_manyAsTopTADs",
                                                            
                                                            "nTADs", "nTopTADs", "nGenes"
                                                            ) ]

all_vars_toplot <- c(all_vars_toplot, "nTop")

# var_to_plot=all_vars_toplot[10]
# cmp_var=all_cmp_var[10]


# for(var_to_plot in all_vars_toplot){
#   if(var_to_plot == "nTop") {
#     cat("A\n")
#     # plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, ".+enes_"), colnames(all_ds_DT))]    
#     plotDT <- all_ds_DT[, grepl("^nTopTADs_genes$|^nTopGenes_manyAsTopTADs$|^nTopGenes_manyAsPercent$", colnames(all_ds_DT))]
#   } else {
#     # cat("B\n")
#     plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, "$"), colnames(all_ds_DT))]
#   }
#   stopifnot(ncol(plotDT) == 3)
#   plotDT$dataset <- all_ds_DT$dataset
#   plotDT_m <- melt(plotDT, id="dataset")
#   plotDT_m$dataset <- factor(plotDT_m$dataset, levels = all_ds_DT$dataset[order(all_ds_DT$aucFCC, decreasing=TRUE)])
#   plotDT_m$variable_leg <- gsub(paste0("_", var_to_plot), "", plotDT_m$variable)
#   plotTit <- paste0(var_to_plot)
#   mySub <- paste0(myGenSub)
# }

for(var_to_plot in all_vars_toplot){
  
  if(var_to_plot == "nTop") {
    cat("A\n")
    # plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, ".+enes_"), colnames(all_ds_DT))]    
    plotDT <- all_ds_DT[, grepl("^nTopTADs_genes$|^nTopGenes_manyAsTopTADs$|^nTopGenes_manyAsPercent$", colnames(all_ds_DT))]
  } else {
    # cat("B\n")
    plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, "$"), colnames(all_ds_DT))]
  }

  plotDT$dataset <- all_ds_DT$dataset
  plotDT_m <- melt(plotDT, id="dataset")
  
  plotDT_m$dataset <- factor(plotDT_m$dataset, levels = all_ds_DT$dataset[order(all_ds_DT$aucFCC, decreasing=TRUE)])
  
  plotDT_m$variable_leg <- gsub(paste0("_", var_to_plot), "", plotDT_m$variable)
  
  # plotTit <- paste0(var_to_plot)
  plotTit <- paste0(GO_aliases_common_top[var_to_plot])
  mySub <- paste0(myGenSub)
  
  # stopifnot(length(unique(plotDT_m$variable)) == 2)
  stopifnot(length(unique(plotDT_m$variable)) == 3)
    
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2", "chartreuse4"), unique(plotDT_m$variable)),
                        # values = setNames(c("dodgerblue4", "darkorange2"), unique(plotDT_m$variable)),
                        labels = setNames(unique(plotDT_m$variable_leg), unique(plotDT_m$variable)))+
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
    
    outFile <- file.path(outFold, paste0(var_to_plot, "_barplot.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    plotDT_m$variable_leg <- as.character(plotDT_m$variable_leg)
    
#    my_comparisons <- list(unique(plotDT_m$variable_leg))
    
    

x <- combn(unique(plotDT_m$variable_leg),2)
my_comparisons <- lapply(seq_len(ncol(x)), function(i) unlist(x[,i]))

myylab <- GO_aliases_common_top[var_to_plot]


    p_box <- ggviolin(plotDT_m, x = "variable_leg", y = "value", 
                      fill = "variable_leg",
                      palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                      title = plotTit,
                      xlab="",
                      ylab = myylab,
                      # sub=mySub,
                      legend.title = "",
                      add = "boxplot", add.params = list(fill = "white"))+
      stat_compare_means(comparisons = my_comparisons, 
                         method="wilcox",
                         # position = "bottomleft",
                         label = "p.format")    + # Add significance levels
      # stat_compare_means()
    
      # stat_compare_means(comparisons = my_comparisons, 
      #                    # position = "bottomleft",
      #                    label = c("p.format"))+ # Add significance levels
      # stat_compare_means(comparisons = my_comparisons,
      #                    # position="bottom",
      #                    label.x.npc = 0,
      #                    label.y.npc=0,
      #                    method="wilcox",
      #                    label = c("p.format"))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot.", plotType))
    ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
}

######################################################################################
######################################################################################
######################################################################################



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

