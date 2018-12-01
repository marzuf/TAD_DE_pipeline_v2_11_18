startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R\n"))

# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 

# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 5 avg
# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 10 avg
# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 15 avg
# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 20 avg
# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 25 avg
# Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 66 avg

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
suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


source("ds_list_run1.R")


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

outFold <- file.path("CMP_TADs_GENES_GO_pvalSelect_top", enricher_ontologyType,nTopDS,rankVar)
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
    
    # run GO analysis for these 3 sets of entrez genes: selectTADs_genes, selectGenes
    #***** 1) selectTADs_genes
    if(nSelectTADs_genes > 0) {
      stopifnot(typeof(universe_tadGenes) == typeof(selectTADs_genes))
      stopifnot(typeof(universe_tadGenes) == typeof(c5_msigdb$gene))
      stopifnot(length(universe_tadGenes) >= length(selectTADs_genes))
      selectTADs_genes_enrich <- enricher(gene = selectTADs_genes, 
                                          TERM2GENE=c5_msigdb,
                                          universe = universe_tadGenes,
                                          pvalueCutoff = enricher_pvalueCutoff, 
                                          pAdjustMethod = enricher_pAdjustMethod, 
                                          minGSSize = enricher_minGSSize, 
                                          maxGSSize = enricher_maxGSSize, 
                                          qvalueCutoff =enricher_qvalueCutoff)
      
      selectTADs_genes_resultDT <- selectTADs_genes_enrich@result
      selectTADs_genes_resultDT <- selectTADs_genes_resultDT[order(selectTADs_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]

      if(nrow(selectTADs_genes_resultDT) < nTopResultsToSave) {
        selectTADs_genes_enrichSaveDT <- selectTADs_genes_resultDT
     } else {
        selectTADs_genes_enrichSaveDT <- selectTADs_genes_resultDT[1:nTopResultsToSave,]
      }
      outFile <- file.path(outFold, paste0(curr_ds,"_selectTADs_genes_enrichSaveDT.Rdata"))
      save(selectTADs_genes_enrichSaveDT, file = outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      
      tmpDT <- selectTADs_genes_resultDT[selectTADs_genes_resultDT[,paste0(padjVarGO)] <= pvalSelectGO,]
      
      if(nrow(tmpDT) > 0) {
        stopifnot(rownames(tmpDT) %in% names(nGenesByGo))
        nGenesGO <- sapply(rownames(tmpDT), function(x) as.numeric(nGenesByGo[x]))
        stopifnot(is.numeric(nGenesGO))
        stopifnot(!is.na(nGenesGO))
        GOdim <- tmpDT$Count/nGenesGO
        stopifnot(GOdim > 0 & GOdim <= 1)
        selectTADs_genes_meanGOdim <- mean(GOdim)
        stopifnot(!is.na(selectTADs_genes_meanGOdim))  
        
      }else {
        selectTADs_genes_meanGOdim <- NA
      }
      

      selectTADs_genes_signifGOterm <- as.character(selectTADs_genes_resultDT$ID[selectTADs_genes_resultDT[,paste0(padjVarGO)] <= pvalSelectGO])
      nSelectTADs_genes_signifGOterm <- length(selectTADs_genes_signifGOterm)
      
      selectTADs_genes_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", selectTADs_genes_signifGOterm)))
      
      length(selectTADs_genes_signifGOterm_tmp)
      # 58
      sum(selectTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
      # 46
      sum(selectTADs_genes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
      # 46
      sum(selectTADs_genes_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
      # 46
      
      selectTADs_genes_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% selectTADs_genes_signifGOterm_tmp]
      nSelectTADs_genes_signifGOid <- length(selectTADs_genes_signifGOid)
      
      # if(nSelectTADs_genes_signifGOterm > 0)
      #   stopifnot(nSelectTADs_genes_signifGOid > 0) # not always TRUE (task 3)
      
    } else {
      selectTADs_genes_resultDT <- data.frame(p.adjust = NA)
      # selectTADs_genes_topGo <- NA
      # selectTADs_genes_nTop <- NA
      selectTADs_genes_signifGOterm <- NA
      selectTADs_genes_signifGOid <- NA
      nSelectTADs_genes_signifGOid <- NA
      nSelectTADs_genes_signifGOterm <- NA
      selectTADs_genes_meanGOdim <- NA
    }
    
    #***** 2) selectGenes
    if(nSelectGenes > 0) {
      stopifnot(typeof(universe_rnaGenes) == typeof(selectGenes))
      stopifnot(typeof(universe_rnaGenes) == typeof(c5_msigdb$gene))
      stopifnot(length(universe_rnaGenes) >= length(selectGenes))
      selectGenes_enrich <- enricher(gene = selectGenes, 
                                     TERM2GENE=c5_msigdb,
                                     universe = universe_rnaGenes,
                                     pvalueCutoff = enricher_pvalueCutoff, 
                                     pAdjustMethod = enricher_pAdjustMethod, 
                                     minGSSize = enricher_minGSSize, 
                                     maxGSSize = enricher_maxGSSize, 
                                     qvalueCutoff =enricher_qvalueCutoff)
      
      selectGenes_resultDT <- selectGenes_enrich@result
      selectGenes_resultDT <- selectGenes_resultDT[order(selectGenes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]


      if(nrow(selectGenes_resultDT) < nTopResultsToSave) {
        selectGenes_enrichSaveDT <- selectGenes_resultDT
     } else {
        selectGenes_enrichSaveDT <- selectGenes_resultDT[1:nTopResultsToSave,]
      }
      outFile <- file.path(outFold, paste0(curr_ds, "_selectGenes_enrichSaveDT.Rdata"))
      save(selectGenes_enrichSaveDT, file = outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      tmpDT <- selectGenes_resultDT[selectGenes_resultDT[,paste0(padjVarGO)] <= pvalSelectGO,]
      
      if(nrow(tmpDT) > 0) {
        stopifnot(rownames(tmpDT) %in% names(nGenesByGo))
        nGenesGO <- sapply(rownames(tmpDT), function(x) as.numeric(nGenesByGo[x]))
        stopifnot(is.numeric(nGenesGO))
        stopifnot(!is.na(nGenesGO))
        GOdim <- tmpDT$Count/nGenesGO
        stopifnot(GOdim > 0 & GOdim <= 1)
        selectGenes_meanGOdim <- mean(GOdim)
        stopifnot(!is.na(selectGenes_meanGOdim))  
      } else {
        selectGenes_meanGOdim <- NA
      }

      
      

      selectGenes_signifGOterm <- as.character(selectGenes_resultDT$ID[selectGenes_resultDT[,paste0(padjVarGO)] <= pvalSelectGO])
      nSelectGenes_signifGOterm <- length(selectGenes_signifGOterm)
      
      selectGenes_signifGOterm_tmp <- tolower(gsub("_", " ", gsub("GO_", "", selectGenes_signifGOterm)))
      
      length(selectGenes_signifGOterm_tmp)
      # 5
      sum(selectGenes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_go$TERM))
      # 5
      sum(selectGenes_signifGOterm_tmp %in% tolower(GOid_GOterm_DT_GOdb$TERM))
      # 5
      sum(selectGenes_signifGOterm_tmp %in% tolower(names(GOterm_GOid_cP)))
      # 5
      
      selectGenes_signifGOid  <- GOterm_GOid_cP[tolower(names(GOterm_GOid_cP)) %in% selectGenes_signifGOterm_tmp]
      nSelectGenes_signifGOid <- length(selectGenes_signifGOid)
      
    } else {
      
      selectGenes_resultDT <- data.frame(p.adjust = NA)
      # selectGenes_topGo <- NA
      # selectGenes_nTop <- NA
      selectGenes_signifGOid <- NA
      selectGenes_signifGOterm <- NA
      nSelectGenes_signifGOid <- NA
      nSelectGenes_signifGOterm <- NA
      selectGenes_meanGOdim <- NA
    }


###### GO analysis

if(!is.na(nSelectTADs_genes_signifGOid) & nSelectTADs_genes_signifGOid > 0) {
  
  selectTADs_genes_signifGOmin <- minimal_set(go, selectTADs_genes_signifGOid)
  nSelectTADs_genes_signifGOmin <- length(selectTADs_genes_signifGOmin)
  selectTADs_genes_gNEL <- inducedGraph(dag = GO_g, startNodes = selectTADs_genes_signifGOmin[selectTADs_genes_signifGOmin %in% nodes(GO_g)])
  
  # !!! need to reverse the graph !!!
  selectTADs_genes_gNEL_rev <- reverseArch(selectTADs_genes_gNEL) # topGO
  selectTADs_genes_gNEL_rev_noAll <- removeNode("all", selectTADs_genes_gNEL_rev) # not needed ? -> do not want this root
  selectTADs_genes_g <- igraph.from.graphNEL(graphNEL = selectTADs_genes_gNEL_rev_noAll)
  
  stopifnot(is.connected(selectTADs_genes_g))
  stopifnot(is.directed(selectTADs_genes_g))
  stopifnot(is.dag(selectTADs_genes_g))
  
  selectTADs_genes_g_nNodes <- length(V(selectTADs_genes_g))
  
  # density
  selectTADs_genes_g_density <- edge_density(selectTADs_genes_g)
  
  # diameter
  selectTADs_genes_g_diameter <- diameter(selectTADs_genes_g)
  
  # mean geodesic distance
  selectTADs_genes_g_meanDist <- mean_distance(selectTADs_genes_g, directed = FALSE) # directed TRUE or FALSE ???
  selectTADs_genes_g_meanDistDir <- mean_distance(selectTADs_genes_g, directed = TRUE) # directed TRUE or FALSE ???

  # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
  selectTADs_genes_g_meanEccentricity <- mean(eccentricity(selectTADs_genes_g))
  
  # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
  selectTADs_genes_g_meanBetweenness <- mean(betweenness(selectTADs_genes_g))    
  
} else {
  selectTADs_genes_signifGOmin <- NA
  nSelectTADs_genes_signifGOmin <- NA
  selectTADs_genes_g_nNodes <- NA
  selectTADs_genes_g_density <- NA
  selectTADs_genes_g_diameter <- NA
  selectTADs_genes_g_meanDist <- NA
  selectTADs_genes_g_meanDistDir <- NA
  selectTADs_genes_g_meanEccentricity <- NA
  selectTADs_genes_g_meanBetweenness <- NA
}
    
    cat("nSelectGenes_signifGOid = ", nSelectGenes_signifGOid, "\n")

if(!is.na(nSelectGenes_signifGOid) & nSelectGenes_signifGOid > 0) {
  
  selectGenes_signifGOmin <- minimal_set(go, selectGenes_signifGOid)
  nSelectGenes_signifGOmin <- length(selectGenes_signifGOmin)
  selectGenes_gNEL <- inducedGraph(dag = GO_g, startNodes = selectGenes_signifGOmin[selectGenes_signifGOmin %in% nodes(GO_g)])
  
  # !!! need to reverse the graph !!!
  selectGenes_gNEL_rev <- reverseArch(selectGenes_gNEL) # topGO
  selectGenes_gNEL_rev_noAll <- removeNode("all", selectGenes_gNEL_rev) # not needed ? -> do not want this root
  selectGenes_g <- igraph.from.graphNEL(graphNEL = selectGenes_gNEL_rev_noAll)
  
  stopifnot(is.connected(selectGenes_g))
  stopifnot(is.directed(selectGenes_g))
  stopifnot(is.dag(selectGenes_g))
  
  selectGenes_g_nNodes <- length(V(selectGenes_g))
  
  # density
  selectGenes_g_density <- edge_density(selectGenes_g)
  
  # diameter
  selectGenes_g_diameter <- diameter(selectGenes_g)
  
  # mean geodesic distance
  selectGenes_g_meanDist <- mean_distance(selectGenes_g, directed = FALSE) # directed TRUE or FALSE ???
  selectGenes_g_meanDistDir <- mean_distance(selectGenes_g, directed = TRUE) # directed TRUE or FALSE ???
  
  # eccentricity - shortest path distance from the farthest other node in the graph. - not relevant ?
  selectGenes_g_meanEccentricity <- mean(eccentricity(selectGenes_g))
  
  # centrality - how many steps is required to access every other vertex from a given vertex - not relevant ?  
  selectGenes_g_meanBetweenness <- mean(betweenness(selectGenes_g))   
    
} else {
  selectGenes_signifGOmin <- NA
  nSelectGenes_signifGOmin <- NA
  selectGenes_g_nNodes <- NA
  selectGenes_g_density <- NA
  selectGenes_g_diameter <- NA
  selectGenes_g_meanDist <- NA
  selectGenes_g_meanDistDir <- NA
  selectGenes_g_meanEccentricity <- NA
  selectGenes_g_meanBetweenness <- NA
  
}
    cat("nSelectGenes_signifGOmin = ", nSelectGenes_signifGOmin, "\n")
    
    nIntersectGenes <- length(intersect(na.omit(selectGenes), na.omit(selectTADs_genes)))
    nUnionGenes <- length(union(na.omit(selectGenes), na.omit(selectTADs_genes)))
  
    nIntersectSignifGOid <- length(intersect(na.omit(selectTADs_genes_signifGOid), na.omit(selectGenes_signifGOid)))
    nIntersectSignifGOterm <- length(intersect(na.omit(selectTADs_genes_signifGOterm), na.omit(selectGenes_signifGOterm)))
    nIntersectSignifGOmin <- length(intersect(na.omit(selectTADs_genes_signifGOmin), na.omit(selectGenes_signifGOmin)))
    
    nSelectTADs_genes_signifGOidRatio <- nSelectTADs_genes_signifGOid/nSelectTADs_genes
    nSelectTADs_genes_signifGOtermRatio <- nSelectTADs_genes_signifGOterm/nSelectTADs_genes
    nSelectTADs_genes_signifGOminRatio <- nSelectTADs_genes_signifGOmin/nSelectTADs_genes
    
    nSelectGenes_signifGOidRatio <- nSelectGenes_signifGOid/nSelectGenes
    nSelectGenes_signifGOtermRatio <- nSelectGenes_signifGOterm/nSelectGenes
    nSelectGenes_signifGOminRatio <- nSelectGenes_signifGOmin/nSelectGenes
    
    intersectSignifGOidRatio <- length(intersect(na.omit(selectTADs_genes_signifGOid), na.omit(selectGenes_signifGOid)))/
      length(union(na.omit(selectTADs_genes_signifGOid), na.omit(selectGenes_signifGOid)))
    intersectSignifGOtermRatio <- length(intersect(na.omit(selectTADs_genes_signifGOterm), na.omit(selectGenes_signifGOterm)))/
      length(union(na.omit(selectTADs_genes_signifGOterm), na.omit(selectGenes_signifGOterm)))
    intersectSignifGOminRatio <- length(intersect(na.omit(selectTADs_genes_signifGOmin), na.omit(selectGenes_signifGOmin)))/
      length(union(na.omit(selectTADs_genes_signifGOmin), na.omit(selectGenes_signifGOmin)))
    
    stopifnot(na.omit(intersectSignifGOidRatio) >= 0 & na.omit(intersectSignifGOidRatio <=1 ))
    stopifnot(na.omit(intersectSignifGOtermRatio) >= 0 & na.omit(intersectSignifGOtermRatio <=1 ))
    stopifnot(na.omit(intersectSignifGOminRatio) >= 0 & na.omit(intersectSignifGOminRatio <=1 ))
    
    selectGenes_intersectRatio <- nIntersectGenes/nSelectGenes
    selectTADs_genes_intersectRatio <- nIntersectGenes/nSelectTADs_genes
    
    intersectGenesRatio <- nIntersectGenes/nUnionGenes
    
    stopifnot(na.omit(selectGenes_intersectRatio) >= 0 & na.omit(selectGenes_intersectRatio <=1 ))
    stopifnot(na.omit(selectTADs_genes_intersectRatio) >= 0 & na.omit(selectTADs_genes_intersectRatio <=1 ))
    stopifnot(na.omit(intersectGenesRatio) >= 0 & na.omit(intersectGenesRatio <=1 ))

    data.frame(
      dataset = curr_ds, 
      
    nTADs = nTADs,
    nSelectTADs = nSelectTADs,
    
    nSelectTADs_genes=nSelectTADs_genes,
    nSelectTADs_genes_signifGOid=nSelectTADs_genes_signifGOid,
    nSelectTADs_genes_signifGOterm=nSelectTADs_genes_signifGOterm,
    nSelectTADs_genes_signifGOmin=nSelectTADs_genes_signifGOmin,
    
    
    selectTADs_genes_g_density=selectTADs_genes_g_density,
    selectTADs_genes_g_diameter=selectTADs_genes_g_diameter,
    selectTADs_genes_g_meanDist=selectTADs_genes_g_meanDist,
    selectTADs_genes_g_meanDistDir=selectTADs_genes_g_meanDistDir,
    selectTADs_genes_g_meanEccentricity=selectTADs_genes_g_meanEccentricity,
    selectTADs_genes_g_meanBetweenness=selectTADs_genes_g_meanBetweenness,
    
    nGenes=nGenes,
    nSelectGenes=nSelectGenes,
    nSelectGenes_signifGOid=nSelectGenes_signifGOid,
    nSelectGenes_signifGOterm=nSelectGenes_signifGOterm,
    nSelectGenes_signifGOmin=nSelectGenes_signifGOmin,
    
    selectGenes_g_density=selectGenes_g_density,
    selectGenes_g_diameter=selectGenes_g_diameter,
    selectGenes_g_meanDist=selectGenes_g_meanDist,
    selectGenes_g_meanDistDir=selectGenes_g_meanDistDir,
    selectGenes_g_meanEccentricity=selectGenes_g_meanEccentricity,
    selectGenes_g_meanBetweenness=selectGenes_g_meanBetweenness,
    
    nIntersectGenes = nIntersectGenes,
    nUnionGenes = nUnionGenes,
    
    nIntersectSignifGOid = nIntersectSignifGOid,
    nIntersectSignifGOterm = nIntersectSignifGOterm,
    nIntersectSignifGOmin = nIntersectSignifGOmin,
    
    nSelectTADs_genes_signifGOidRatio = nSelectTADs_genes_signifGOidRatio,
    nSelectTADs_genes_signifGOtermRatio = nSelectTADs_genes_signifGOtermRatio,
    nSelectTADs_genes_signifGOminRatio = nSelectTADs_genes_signifGOminRatio,
    
    nSelectGenes_signifGOidRatio = nSelectGenes_signifGOidRatio,
    nSelectGenes_signifGOtermRatio = nSelectGenes_signifGOtermRatio,
    nSelectGenes_signifGOminRatio = nSelectGenes_signifGOminRatio,
    
    intersectSignifGOidRatio = intersectSignifGOidRatio,
    intersectSignifGOtermRatio = intersectSignifGOtermRatio,
    intersectSignifGOminRatio = intersectSignifGOminRatio,
    
    
    selectGenes_intersectRatio = selectGenes_intersectRatio,
    selectTADs_genes_intersectRatio = selectTADs_genes_intersectRatio,
    
    intersectGenesRatio = intersectGenesRatio,
    
    selectTADs_genes_meanGOdim = selectTADs_genes_meanGOdim,
    selectGenes_meanGOdim = selectGenes_meanGOdim,
    
    stringsAsFactors=FALSE
    
    )

    #save(enricher_all_results, file = outFile)
    #cat(paste0("... written: ", outFile, "\n"))
    
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

stopifnot(all_ds_DT$dataset %in% names(all_ds_aucFCC))
all_ds_DT$aucFCC <- all_ds_aucFCC[all_ds_DT$dataset]

stopifnot(all_ds_DT$dataset %in% names(all_ds_aucCoexprDist))
all_ds_DT$aucCoexprDist <- all_ds_aucCoexprDist[all_ds_DT$dataset]


stopifnot(!duplicated(all_ds_DT$dataset))

stopifnot(na.omit(all_ds_DT$intersectSignifGOidRatio) >= 0 & na.omit(all_ds_DT$intersectSignifGOidRatio <=1 ))
stopifnot(na.omit(all_ds_DT$intersectSignifGOtermRatio) >= 0 & na.omit(all_ds_DT$intersectSignifGOtermRatio <=1 ))
stopifnot(na.omit(all_ds_DT$intersectSignifGOminRatio) >= 0 & na.omit(all_ds_DT$intersectSignifGOminRatio <=1 ))

stopifnot(na.omit(all_ds_DT$selectGenes_intersectRatio) >= 0 & na.omit(all_ds_DT$selectGenes_intersectRatio <=1 ))
stopifnot(na.omit(all_ds_DT$selectTADs_genes_intersectRatio) >= 0 & na.omit(all_ds_DT$selectTADs_genes_intersectRatio <=1 ))
stopifnot(na.omit(all_ds_DT$intersectGenesRatio) >= 0 & na.omit(all_ds_DT$intersectGenesRatio <=1 ))


##########################################################################################
# save table in text file
outDT <- all_ds_DT
outDT <- outDT[order(outDT$aucFCC, decreasing=TRUE),]
write.table(outDT, file = file.path(outFold, "all_ds_DT.txt"), sep="\t", quote=F, col.names=T, row.names=F)
# stop("..ok..\n")
##########################################################################################
# set title name for the different variables
stopifnot(colnames(all_ds_DT) %in% names(GO_aliases_pvalSelect)) # loaded from analysis_utils.R

# subtitle for the barplots:
myGenSub <- paste0("pvalSelect = ", pvalSelect, "; pvalSelectGO = ", pvalSelectGO)


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
  all_VDs <- list()
  
  ds_dt <- all_ds_DT[all_ds_DT$dataset == curr_ds,]
  stopifnot(nrow(ds_dt) == 1)
  
  all_VDs[["nGenes"]] <- draw_myVenn(var1="nSelectGenes", 
              var2="nSelectTADs_genes", 
              varIntersect = "nIntersectGenes",
              tit = curr_ds,
              subTit = paste0(GO_aliases_common_pvalSelect[gsub("Genes", "", "nSelectGenes")], " (nSelectTADs = ", ds_dt$nSelectTADs, ")"),
              dataDT = ds_dt)
  
  curr_var="signifGOid"
  for(curr_var in c("signifGOid", "signifGOterm", "signifGOmin")) {
    all_VDs[[paste0(curr_var)]] <- draw_myVenn(var1=paste0("nSelectGenes_", curr_var),
                var2=paste0("nSelectTADs_genes_", curr_var),
                varIntersect = paste0("nIntersect", capitalize(curr_var)),
                tit = curr_ds,
                subTit = GO_aliases_common_pvalSelect[paste0(curr_var)],
                dataDT = ds_dt,
                subPatt1 = paste0("_", curr_var),
                subPatt2 = paste0("_", curr_var)
                )
  } # end iterating over variables
  titGrob <- textGrob(paste0(curr_ds),gp=gpar(fontsize=20, fontface="bold"))
  all_vds <- do.call(arrangeGrob, c(all_VDs, nrow=2))
  finalVD <- arrangeGrob(all_vds, ncol = 1, top = titGrob)
  
  outFile <- file.path(outFold, paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_intersect_", curr_ds, "_venn_diagram.", plotType))
  ggsave(finalVD, file = outFile, height = vdHeight, width = vdWidth)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
} # end iterating over datasets


#cat("HELLO4\n")

# plot all datasets

######################################################################################
###################################################################################### BARPLOT INTERSECT
######################################################################################

intersect_dt <- all_ds_DT[, c(which(colnames(all_ds_DT) == "dataset"), grep("intersect.+Ratio", colnames(all_ds_DT)))]
intersect_dt_m <- melt(intersect_dt, id="dataset")
stopifnot(intersect_dt_m$variable %in% names(GO_aliases_pvalSelect))
intersect_dt_m$variable_leg <- GO_aliases_pvalSelect[as.character(intersect_dt_m$variable)]
intersect_dt_m$variable_leg <- gsub(" - intersect", "", intersect_dt_m$variable_leg)

plotTit_intersect <- "Intersect ratio"
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

outFile <- file.path(outFold, paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_all_ds_intersect", "_barplot.", plotType))
ggsave(plot = p_intersect, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
###################################################################################### BARPLOT P VALUES - 1 TYPE SEPARATELY
######################################################################################

# SEE ABOVE FOR SETTING: barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
# barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)

# CMP_TADs_GENES_GO_top/20/TCGAluad_luad_mutKRAS_topGenes_manyAsPercent_enrichSaveDT.Rdata
curr_ds="TCGAcrc_msi_mss"
curr_type="selectGenes"

for(curr_ds in unique(all_ds_DT$dataset)) {
  for(curr_type in c("selectGenes", "selectTADs_genes")) {
    
  resultDT_file <- file.path(outFold, paste0(curr_ds, "_",curr_type, "_enrichSaveDT.Rdata"))
    
    cat(paste0("resultDT_file = ", resultDT_file, "\n"))
    
    if(!file.exists(resultDT_file)) next
    
    stopifnot(file.exists(resultDT_file))
    resultDT <- eval(parse(text = load(resultDT_file)))    
    
    stopifnot(nrow(resultDT) == nTopResultsToSave)
    
    # resultDT$log10_pval <- -log10(resultDT$qvalue)
    resultDT$log10_pval <- -log10(resultDT[,paste0(padjVarGO)])
    # resultDT$geneRatio <- eval(parse(text = resultDT$GeneRatio))
    resultDT$geneRatio <- as.numeric(sapply(resultDT$GeneRatio, function(x) {gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    # resultDT$bgRatio <- eval(parse(text = resultDT$BgRatio))
    resultDT$bgRatio <- as.numeric(sapply(resultDT$BgRatio, function(x) {gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    stopifnot(is.numeric(resultDT$geneRatio))
    stopifnot(is.numeric(resultDT$bgRatio))
    resultDT$foldEnrichment <- resultDT$geneRatio/resultDT$bgRatio
    
    cat("hello1\n")
    
    for(plot_var in barplot_vars) {
      resultDT <- resultDT[order(resultDT[, plot_var], decreasing=TRUE),]
      
      cat("plot_var=", plot_var, "\n")
      cat("colnames=", paste0(colnames(resultDT), collapse=";"), "\n")
      
      stopifnot(plot_var %in% colnames(resultDT))
      myTit <- paste0(barplot_vars_tit[plot_var], " - ", curr_ds, " - ", curr_type)
      outFile <- file.path(outFold,paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_", curr_ds, "_", curr_type, "_top", nTopResultsToSave, "_GO", plot_var, ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(resultDT[,plot_var], 
              main = myTit, 
              # ylab = paste0("-log10(", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[plot_var]),
              names.arg = gsub("_", " ", gsub("GO_", "", resultDT$ID)), 
              las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
}


######################################################################################
###################################################################################### BARPLOT P VALUES - ALL TOGETHER - ITERATE OVER DATASETS - ADDED
######################################################################################

# CMP_TADs_GENES_GO_top/20/TCGAluad_luad_mutKRAS_topGenes_manyAsPercent_enrichSaveDT.Rdata
curr_ds="TCGAcrc_msi_mss"
curr_type="selectGenes"

barplotSubtit <- paste0("pvalSelect=", pvalSelect, "; pvalSelectGO=", pvalSelectGO)

all_ds_pval_DT <- foreach(curr_ds = unique(all_ds_DT$dataset), .combine='rbind') %dopar% {
  ds_pval_dt <- foreach(curr_type = c("selectGenes", "selectTADs_genes"), .combine='rbind') %do% {
    
    resultDT_file <- file.path(outFold, paste0(curr_ds, "_",curr_type, "_enrichSaveDT.Rdata"))
    cat(paste0("resultDT_file = ", resultDT_file, "\n"))
    
    if(!file.exists(resultDT_file)) {
      # result_DT <- data.frame(
      #   ID = NULL,
      #   # Description = NULL,
      #   # GeneRatio = NULL,
      #   # BgRatio = NULL,
      #   # pvalue = NULL,
      #   # p.adjust = NULL,
      #   # qvalue = NULL,
      #   # geneID = NULL,
      #   # Count = NULL, 
      #   log10_pval = NULL,
      #   rank = NULL,
      #   geneType = NULL,
      #   stringsAsFactors=F
      # )
      result_DT <- data.frame(stringsAsFactors = F)
      return(result_DT)
    }
    
    
    cat("hello1\n")
    
    stopifnot(file.exists(resultDT_file))
    resultDT <- eval(parse(text = load(resultDT_file)))    
    stopifnot(nrow(resultDT) == nTopResultsToSave)
    # resultDT$log10_pval <- -log10(resultDT$qvalue)
    resultDT$log10_pval <- -log10(resultDT[,paste0(padjVarGO)])
    # resultDT$geneRatio <- eval(parse(text = resultDT$GeneRatio))
    resultDT$geneRatio <- as.numeric(sapply(resultDT$GeneRatio, function(x) {gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    # resultDT$bgRatio <- eval(parse(text = resultDT$BgRatio))
    resultDT$bgRatio <- as.numeric(sapply(resultDT$BgRatio, function(x) {gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    stopifnot(is.numeric(resultDT$geneRatio))
    stopifnot(is.numeric(resultDT$bgRatio))
    resultDT$foldEnrichment <- resultDT$geneRatio/resultDT$bgRatio
    
    for(var_plot in barplot_vars) {
      resultDT <- resultDT[order(resultDT[var_plot], decreasing=T),]
      stopifnot(nrow(resultDT) > 0)
      resultDT[,paste0("rank", "_", var_plot)] <- 1:nrow(resultDT) # not use rank because of ties
    }
    
    
              # resultDT <- resultDT[order(resultDT$log10_pval, decreasing=T),]
              # # resultDT <- resultDT[, c("ID", "log10_pval")]
              # resultDT <- resultDT[, c("ID", barplot_vars)]
              # stopifnot(nrow(resultDT) > 0)
              # resultDT$rank <- 1:nrow(resultDT)
  
    # resultDT <- resultDT[, c("ID", barplot_vars, paste0("rank_", barplot_vars))]
    resultDT$geneType <- curr_type
    rownames(resultDT) <- NULL
    resultDT
  } # end iterate selectGenes, selectTADs_genes
  
  
  for(var_plot in barplot_vars) {
    
    ds_pval_dt$rank <- ds_pval_dt[,paste0("rank_", var_plot)]
    ds_pval_dt <- ds_pval_dt[order(ds_pval_dt$rank, ds_pval_dt$geneType),]
    ds_pval_dt$x_ID <- c(1:nrow(ds_pval_dt))
    ds_pval_dt$x_ID <- factor(ds_pval_dt$x_ID, levels = unique(as.character(ds_pval_dt$x_ID)))
    ds_pval_dt$x_labels <- gsub("_", " ", gsub("GO_", "", as.character(ds_pval_dt$ID)))
    
    p_DS <- ggplot(ds_pval_dt, aes_string(x = "x_ID", y = var_plot, fill = "geneType")) +
      facet_grid(~rank, switch="x", scale="free")+
      geom_bar(stat="identity", position="dodge") +
      scale_x_discrete(name="", labels=ds_pval_dt$x_labels, breaks = ds_pval_dt$x_ID)   +
      # scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
    scale_y_continuous(name=paste0(barplot_vars_tit[var_plot]),
                         breaks = scales::pretty_breaks(n = 5))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
      coord_cartesian(expand = FALSE) +
      scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
                        labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
      labs(fill  = "") +
      ggtitle(label = paste0(barplot_vars_tit[var_plot], " - top GO - ", curr_ds), subtitle = barplotSubtit) +
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
    outFile <- file.path(outFold,paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_", curr_ds, "_", curr_type, "_top", nTopResultsToSave, "_GO", var_plot, "_selectGenes_selectTADs.", plotType))
    ggsave(p_DS, filename = outFile, height = myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    ds_pval_dt$rank <- NULL
  }

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


for(var_plot in barplot_vars) {
  
  all_ds_pval_DT$rank <- all_ds_pval_DT[,paste0("rank_", var_plot)]
  
  all_ds_pval_DT$rank <- factor(all_ds_pval_DT$rank, levels = as.character(sort(unique(as.numeric(as.character(all_ds_pval_DT$rank))))))
  
  
  p_all <- ggplot(all_ds_pval_DT, aes_string(x = paste0("rank"), y = var_plot, fill = "geneType", col = "geneType")) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter() +
    scale_x_discrete(name="GO rank")+
    # scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
    scale_y_continuous(name=paste0(barplot_vars_tit[var_plot]),
                       breaks = scales::pretty_breaks(n = 5))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
    # coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
                      labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
    scale_colour_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
                        labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"), guide = F)+
    labs(fill  = "") +
    ggtitle(label = paste0(barplot_vars_tit[var_plot], " - top GO - all DS (n=", length(unique(all_ds_pval_DT$dataset)), ")"), subtitle = barplotSubtit) +
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
  
  outFile <- file.path(outFold,paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_allDS_typeComp", "_top", nTopResultsToSave, "_GO", var_plot, "_boxplot.", plotType))
  ggsave(p_all, filename = outFile, height = myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  p_all <- ggplot(all_ds_pval_DT, aes_string(x = paste0("rank"), y = var_plot, fill = "geneType")) +
    geom_boxplot(outlier.shape=NA) +
    #geom_jitter() +
    scale_x_discrete(name="GO rank")+
    # scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
    scale_y_continuous(name=paste0(barplot_vars_tit[var_plot]),
                       breaks = scales::pretty_breaks(n = 5))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
    # coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
                      labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
    labs(fill  = "") +
    ggtitle(label = paste0(barplot_vars_tit[var_plot], " - top GO - all DS (n=", length(unique(all_ds_pval_DT$dataset)), ")"), subtitle = barplotSubtit) +
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
  
  outFile <- file.path(outFold,paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_allDS_typeComp", "_top", nTopResultsToSave, "_GO", var_plot, "_boxplot_nojitter.", plotType))
  ggsave(p_all, filename = outFile, height = myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  all_ds_pval_DT$rank <- NULL
  
}


######################################################################################
###################################################################################### SCATTERPLOT AUC FCC VS. AUC COEXPR DIST.
######################################################################################

# scatterplot vs. aucFCC and aucCoexprDist

          # outFile <- file.path(outFold, "all_ds_DT.Rdata")
          # stopifnot(file.exists(outFile))
          # all_ds_DT <- eval(parse(text = load(outFile)))
          # 
          # stopifnot(all_ds_DT$dataset %in% names(all_ds_aucFCC))
          # all_ds_DT$aucFCC <- all_ds_aucFCC[all_ds_DT$dataset]
          # 
          # stopifnot(all_ds_DT$dataset %in% names(all_ds_aucCoexprDist))
          # all_ds_DT$aucCoexprDist <- all_ds_aucCoexprDist[all_ds_DT$dataset]

stopifnot(as.character(all_ds_DT$dataset)  %in% names(dataset_proc_colors) )

all_vars <- colnames(all_ds_DT)
all_vars <- all_vars[all_vars != "dataset"]

ref_vars <- c("aucFCC", "aucCoexprDist")

stopifnot(all_vars %in% names(GO_offSets_pvalSelect))
stopifnot(ref_vars %in% names(GO_offSets_pvalSelect))
stopifnot(all_vars %in% names(GO_legPos_pvalSelect))
stopifnot(ref_vars %in% names(GO_legPos_pvalSelect))

ref_var=ref_vars[1]
curr_var=all_vars[6]

for(ref_var in ref_vars){
  for(curr_var in all_vars) {
    
    if(ref_var == curr_var) next
    
    if(grepl("selectTADs_genes", curr_var) | grepl("SelectTADs_genes", curr_var)) {
      mysub <- paste0("selectTADs_genes - ", myGenSub)
    }else if(grepl("selectGenes", curr_var) | grepl("SelectGenes", curr_var)) {
      mysub <- paste0("selectGenes - ", myGenSub)
    } else {
      mysub <- paste0(myGenSub)
    }
    
    # myylab <- unique(gsub("selectTADs_genes_", "", 
    #                         gsub("selectGenes_", "", 
    #                              gsub("nSelectGenes_", "", 
    #                                   gsub("nSelectTADs_genes_", "",  curr_var)))))
    # 
    myylab <- GO_aliases_pvalSelect[curr_var]
    myxlab <- GO_aliases_pvalSelect[ref_var]
    
    mymain <- paste0(myylab, " vs. ", ref_var)
    
    myx <- all_ds_DT[,ref_var]
    myy <- all_ds_DT[,curr_var]
    
    if(all(is.na(myx)) | all(is.na(myy))) next
    
    if( sum(!is.na(myy)) < 2 | sum(!is.na(myx)) < 2 ) next # need at least 2 to draw the line
    
    curr_colors <- dataset_proc_colors[as.character(all_ds_DT$dataset)]
    
    outFile <- file.path(outFold, paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_", curr_var, "_" , ref_var, ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(x=myx,
         y=myy,
         main = paste0(myylab, " vs. ", ref_var),
         # xlab =  paste0(ref_var),
         xlim = range(myx, na.rm=T) + c(-GO_offSets_pvalSelect[ref_var], GO_offSets_pvalSelect[ref_var]),
         ylim = range(myy, na.rm=T) + c(-GO_offSets_pvalSelect[curr_var], GO_offSets_pvalSelect[curr_var]),
         xlab = myxlab,
         ylab = myylab,
         pch=16,cex=0.7)
    add_curv_fit(x=myx,
                 y=myy, withR2 = F, lty=2)
    addCorr(x=myx,
            y=myy, legPos = "topleft", corMet="spearman", bty="n")
    text(x = myx, 
         y = myy, 
         labels = all_ds_DT[,"dataset"], 
         col = curr_colors,
         cex=0.7)
    
    legend(
           GO_legPos_pvalSelect[curr_var], 
           #"bottomright",
           legend=names(my_colors),
           lty=-1,
           pch=16,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.6)
    
    mtext(mysub,side=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

################################################################################### 

all_cmp_var <- c("selectGenes", "selectTADs_genes")

all_vars_toplot <- unique(gsub("selectTADs_genes_", "",
                      gsub("selectGenes_", "",
                           gsub("nSelectGenes_", "",
                                gsub("nSelectTADs_genes_", "",  colnames(all_ds_DT))))))
# all_vars_toplot <- colnames(all_ds_DT)
all_vars_toplot <- all_vars_toplot[! all_vars_toplot %in% c("dataset", "aucFCC", "aucCoexprDist",
                                                            "nSelectGenes", "nSelectTADs_genes",
                                                            
                                                            "nIntersectSignifGOterm","nIntersectSignifGOid", "nIntersectSignifGOmin",
                                                            "intersectSignifGOtermRatio","intersectSignifGOidRatio", "intersectSignifGOminRatio",
                                          
                                                            "nIntersectGenes", "nUnionGenes",
                                                            "intersectGenesRatio",
                                                            
                                                            "nTADs", "nSelectTADs", "nGenes"
                                                            ) ]
all_vars_toplot <- c(all_vars_toplot, "nSelect")

var_to_plot=all_vars_toplot[1]
cmp_var=all_cmp_var[1]

all_vars_toplot[ ! all_vars_toplot %in% names(GO_aliases_common_pvalSelect) ]

stopifnot(all_vars_toplot %in% names(GO_aliases_common_pvalSelect))

for(var_to_plot in all_vars_toplot){
  
  # cat("var_to_plot = ", var_to_plot, "\n")
  # cat(paste0(colnames(all_ds_DT), collapse=","), "\n")
  
  if(var_to_plot == "nSelect") {
    plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, ".+enes$"), colnames(all_ds_DT))]    
  } else {
    plotDT <- all_ds_DT[, grepl(paste0(var_to_plot, "$"), colnames(all_ds_DT))]
  }

  plotDT$dataset <- all_ds_DT$dataset
  plotDT_m <- melt(plotDT, id="dataset")
  
  plotDT_m$dataset <- factor(plotDT_m$dataset, levels = all_ds_DT$dataset[order(all_ds_DT$aucFCC, decreasing=TRUE)])
  
  plotDT_m$variable_leg <- gsub(paste0("_", var_to_plot), "", plotDT_m$variable)
  # stopifnot(plotDT_m$variable_leg %in% names(GO_aliases_common_pvalSelect))
  # plotDT_m$variable_leg <- GO_aliases_common_pvalSelect[as.character(plotDT_m$variable_leg)]
  
  # plotTit <- paste0(var_to_plot)
  plotTit <- paste0(GO_aliases_common_pvalSelect[var_to_plot])
  mySub <- paste0(myGenSub)
    
  stopifnot(length(unique(plotDT_m$variable)) == 2)
    
    p_var <- ggplot(plotDT_m, aes(x = dataset, y = value, fill = variable_leg)) +
      ggtitle(plotTit, subtitle = mySub)+
      geom_bar(stat="identity", position = "dodge")+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(""),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_fill_manual(values = setNames(c("dodgerblue4", "darkorange2"), unique(plotDT_m$variable_leg)),
                        labels = setNames(unique(plotDT_m$variable_leg), unique(plotDT_m$variable_leg)))+
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
    outFile <- file.path(outFold, paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_", var_to_plot, "_barplot.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    plotDT_m$variable_leg <- as.character(plotDT_m$variable_leg)
    
    # my_comparisons <- list(c(ref_var, cmp_var))
    my_comparisons <- list(unique(plotDT_m$variable_leg))
    # my_comparisons <- list(unique(plotDT_m$variable_leg))
    # my_comparisons <- list(c("nSelectTADs_genes", "nSelectGenes"))
    # cat(my_comparisons)
    # myylab <- paste0(var_to_plot)
    myylab <- GO_aliases_common_pvalSelect[var_to_plot]
    # myylab <- paste0(GO_aliases_pvalSelect[var_to_plot])
    
    mysubtit <- paste0("all datasets (n=", length(unique(plotDT_m$dataset)), ")")
    
    p_box <- ggviolin(plotDT_m, x = "variable_leg", y = "value", 
                      fill = "variable_leg",
                      palette = c("#00AFBB", "#FC4E07"),
                      title = plotTit,
                      subtitle = mysubtit,
                      xlab="",
                      ylab = myylab,# gsub("_", " ", var_to_plot),
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
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

   
    outFile <- file.path(outFold, paste0("GOtype", enricher_ontologyType, "_nTopDS", nTopDS, "_rankVar", rankVar, "_", var_to_plot, "_violinplot.", plotType))
    ggsave(plot = p_box, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
}


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

