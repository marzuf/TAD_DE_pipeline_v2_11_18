

######################################################################################################################################################################################################
###################################################################################################################################################################################################### add_curv_fit (function)
######################################################################################################################################################################################################

add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### addCorr (function)
######################################################################################################################################################################################################

addCorr <- function(x, y, legPos="topright", corMet="pearson", ...) {
  corMet <- tolower(corMet)
  stopifnot(corMet %in% c("pearson", "kendall", "spearman"))
  x2 <- x[!is.na(x) & !is.na(y)]
  y2 <- y[!is.na(x) & !is.na(y)]
  x <- x2
  y <- y2
  stopifnot(length(x) == length(y))

  if(length(x) < 3) {
    legTxt <- paste0(paste0(toupper(substr(corMet,1,1)), "CC"), " = NA", "\n", "(# obs. < 3)")
  } else {
    ct <- cor.test(x,y, method = corMet)
    corCoeff <- ct$estimate
    corPval <- ct$p.value
    legTxt <- paste0(paste0(toupper(substr(corMet,1,1)), "CC"), " = ", round(corCoeff, 4), "\n", "(p-val = ", sprintf("%2.2e", corPval), ")")
  }
  legend(legPos, legend = legTxt, ...)
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### printAndLog (function)
######################################################################################################################################################################################################

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### densplot (function)
######################################################################################################################################################################################################

densplot <- function(x,y, pch=19, cex=1, ...){
	df <- data.frame(x,y)
	d <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
	df$dens <- col2rgb(d)[1,] + 1L
	cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$dens]
	df <- df[order(df$dens),]
	plot(df$x,df$y, pch=pch, col=df$col, ...)
}


######################################################################################################################################################################################################
###################################################################################################################################################################################################### plot_multiDens (function)
######################################################################################################################################################################################################

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### my_plot_function (function)
######################################################################################################################################################################################################

my_plot_function <- function(varX, varY, mylabels = NULL, withLab=F, ...) {
  myx <- varX
  myy <- varY
  plot(x = myx,
       y = myy,
#       xlab = ifelse(grepl("auc", var1), var1, paste0("signif. GO - ", var1)),
#       ylab = paste0(var2),
#       main = paste0(var2, " vs. ", var1),
       pch=16, cex=0.7,
        ...
  )
  add_curv_fit(x=myx,
               y=myy, withR2 = F, lty=2)
  addCorr(x=myx,
          y=myy, bty="n")
  if(withLab){
    stopifnot(!is.null(mylabels))
    text(x = myx, y=myy, labels = mylabels, cex=0.6)
  }
}




######################################################################################################################################################################################################
###################################################################################################################################################################################################### plot_cumsumDiff05_returnVect (function)
######################################################################################################################################################################################################

plot_cumsumDiff05_returnVect <- function(observ_vect, permut_DT, pointObsCol = "black", 
           my_stat = "ratioDown",
           polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),
           departureValue = 0.5, drawline=FALSE, ...) {
  
  plotTrueType <- ifelse(drawline, "l", "p")

#  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
#  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
#  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)

  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  
  plot(cumsum(abs(observ_vect - departureValue)) ~ x_val,
#       main= my_main,
#       cex.main = cexMain,
#       xlab= my_xlab, 
#       ylab= my_ylab,
       type = plotTrueType,
       pch = 16, cex = 0.7,
       bty="l", ...)
  polygon(x = c(x_val, rev(x_val)), 
          y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
          border=NA,
          col = polygonPermutCol)
  legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
         pt.cex = c(0.7, 2),
         col = c(pointObsCol, polygonPermutCol),
         bty="n")

  return(cumsum(abs(observ_vect - departureValue)))
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### CANCER SUBTYPES (vector)
######################################################################################################################################################################################################

cancer_subAnnot <- c(
"TCGAstad_msi_gs" = "subtypes",
"GSE79209_dysp_nodysp" = "lesions",
"GSE58135_tripleNeg_adjTripleNeg" = "vs_normal",
"GSE40419_normal_cancer" = "vs_normal",
"TCGAstad_EBVneg_EBVpos" = "subtypes",
"GSE77509_normal_ptt" = "vs_normal",
"GSE77509_normal_tumor" = "vs_normal",
"GSE77314_normal_tumor" = "vs_normal",
"GSE71119_dediffSM_MFSM" = "subtypes",
"GSE58135_ERpos_adjERpos" = "vs_normal",
"GSE58135_ERpos_tripleNeg" = "subtypes",
"GSE71119_undiffSM_LMSM" = "subtypes",
"TCGAlaml_laml_mutFLT3" = "mutation",
"TCGAthca_thca_mutBRAF" = "mutation",
"TCGAlihc_lihc_mutCTNNB1" = "mutation",
"GSE77509_ptt_tumor" = "vs_other",
"GSE87340_ad_nl" = "vs_normal",
"GSE81089_normal_nsclc" = "vs_normal",
"TCGApaad_paad_mutKRAS" = "mutation",
"TCGAucec_msi_cnl" = "subtypes",        
"TCGAbrca_lum_bas" = "subtypes",
"TCGAskcm_skcm_mutBRAF" = "mutation",
"TCGAluad_luad_mutKRAS" = "mutation",
"TCGAskcm_skcm_mutCTNNB1" = "mutation",        
"TCGAacc_acc_mutCTNNB1" = "mutation",
"GSE74927_neg_pos" = "subtypes",
"TCGAcrc_msi_mss" =  "subtypes",
"GSE102073_stic_nostic" = "lesions",

"TCGAcesc_adeno_squam" = "subtypes",
"TCGAhnsc_HPVneg_HPVpos" = "subtypes",
"TCGAlgg_IDHwt_IDHmutnc" = "subtypes",
"TCGAsarc_ddlps_lms" = "subtypes",
"TCGAsarc_ddlps_mfs" = "subtypes",
"TCGAsarc_lms_mfs" = "subtypes",
"TCGAtgct_sem_nonsem" = "subtypes",
"TCGAcoad_msi_mss" = "subtypes",

"TCGAskcm_lowinf_highInf" = "subtypes",
"TCGAluad_nonsmoker_smoker" = "subtypes",

"TCGAblca_norm_blca" = "vs_normal",
"TCGAkich_norm_kich" = "vs_normal",
"TCGAlusc_norm_lusc" = "vs_normal",
"TCGAstad_norm_gs" = "vs_normal",



"TCGAacc_wt_mutCTNNB1" = "mutation",
"TCGAgbm_classical_mesenchymal"="subtypes",
"TCGAgbm_classical_neural"="subtypes",
"TCGAgbm_classical_proneural"  ="subtypes",
"TCGAlaml_wt_mutFLT3" = "mutation",
"TCGAlihc_wt_mutCTNNB1"         = "mutation",
"TCGAluad_mutKRAS_mutEGFR"="subtypes",
"TCGAluad_wt_mutKRAS"           = "mutation",
"TCGApaad_wt_mutKRAS" = "mutation",
"TCGAskcm_lowInf_highInf"    ="subtypes",  
"TCGAskcm_wt_mutBRAF" = "mutation",
"TCGAskcm_wt_mutCTNNB1" = "mutation",        
"TCGAstad_EBVpos_gs"="subtypes",
"TCGAthca_mut.RAS_mutBRAF"   =  "subtypes",
"TCGAthca_wt_mutBRAF"           = "mutation"



)

cancer_subColors <- c(
subtypes = "green4",
lesions = "chocolate1",
vs_normal = "blue3",
vs_other = "orchid1",
mutation = "red"
)

stopifnot(cancer_subAnnot %in% names(cancer_subColors))

# to get the color: cancer_subColors[cancer_subAnnot[curr_ds]]

######################################################################################################################################################################################################
###################################################################################################################################################################################################### MAPPING ALIASES FOR GO VARIABLES (nTop) (vector)
######################################################################################################################################################################################################
GO_aliases_top <- c(
  "dataset" = "dataset",
  "nTADs" = "# TADs",
  "nTopTADs" = "# selectTADs",
  "nTopTADs_genes" = "# topTADs_genes",        
  "nTopTADs_genes_signifGOid" = "# signif. GO (id) - topTADs_genes",
  "nTopTADs_genes_signifGOterm" = "# signif. GO (term) - topTADs_genes",
  "nTopTADs_genes_signifGOmin" = "# signif. GO (min. set) - topTADs_genes",
  "topTADs_genes_g_density" = "graph density - topTADs_genes",
  "topTADs_genes_g_diameter" = "graph diameter - topTADs_genes",
  "topTADs_genes_g_meanDist" = "graph mean dist. (undir.) - topTADs_genes",
  "topTADs_genes_g_meanDistDir" = "graph mean dist. (dir.) - topTADs_genes",
  "topTADs_genes_g_meanEccentricity" = "graph mean eccentricity - topTADs_genes",
  "topTADs_genes_g_meanBetweenness" = "graph mean betweenness - topTADs_genes",
  
  "nGenes" = "# genes",
  "nTopGenes_manyAsPercent" = "# topGenes_manyAsPercent",
  "nTopGenes_manyAsPercent_signifGOid" = "# signif. GO (id) - topGenes_manyAsPercent",
  "nTopGenes_manyAsPercent_signifGOterm" = "# signif. GO (term) - topGenes_manyAsPercent",
  "nTopGenes_manyAsPercent_signifGOmin" = "# signif. GO (min. set) - topGenes_manyAsPercent",
  "topGenes_manyAsPercent_g_density" = "graph density - topGenes_manyAsPercent",
  "topGenes_manyAsPercent_g_diameter"= "graph diameter - topGenes_manyAsPercent",
  "topGenes_manyAsPercent_g_meanDist" = "graph mean dist. (undir.) - topGenes_manyAsPercent",
  "topGenes_manyAsPercent_g_meanDistDir" = "graph mean dist. (dir.) - topGenes_manyAsPercent",
  "topGenes_manyAsPercent_g_meanEccentricity" = "graph mean eccentricity - topGenes_manyAsPercent",
  "topGenes_manyAsPercent_g_meanBetweenness" ="graph mean betweenness - topGenes_manyAsPercent",

  "nTopGenes_manyAsTopTADs" = "# topGenes_manyAsTopTADs",
  "nTopGenes_manyAsTopTADs_signifGOid" = "# signif. GO (id) - topGenes_manyAsTopTADs",
  "nTopGenes_manyAsTopTADs_signifGOterm" = "# signif. GO (term) - topGenes_manyAsTopTADs",
  "nTopGenes_manyAsTopTADs_signifGOmin" = "# signif. GO (min. set) - topGenes_manyAsTopTADs",
  "topGenes_manyAsTopTADs_g_density" = "graph density - topGenes_manyAsTopTADs",
  "topGenes_manyAsTopTADs_g_diameter"= "graph diameter - topGenes_manyAsTopTADs",
  "topGenes_manyAsTopTADs_g_meanDist" = "graph mean dist. (undir.) - topGenes_manyAsTopTADs",
  "topGenes_manyAsTopTADs_g_meanDistDir" = "graph mean dist. (dir.) - topGenes_manyAsTopTADs",
  "topGenes_manyAsTopTADs_g_meanEccentricity" = "graph mean eccentricity - topGenes_manyAsTopTADs",
  "topGenes_manyAsTopTADs_g_meanBetweenness" ="graph mean betweenness - topGenes_manyAsTopTADs",

  "nIntersectSignifGOid_manyAsPercent" = "# signif. GO (id) - intersect manyAsPercent",
  "nIntersectSignifGOterm_manyAsPercent" = "# signif. GO (term) - intersect manyAsPercent",
  "nIntersectSignifGOmin_manyAsPercent" = "# signif. GO (min. set) - intersect manyAsPercent",
  "intersectSignifGOidRatio_manyAsPercent" = "ratio signif. GO (id) - intersect manyAsPercent",
  "intersectSignifGOtermRatio_manyAsPercent" = "ratio signif. GO (term) - intersect manyAsPercent",
  "intersectSignifGOminRatio_manyAsPercent" = "ratio signif. GO (min. set) - intersect manyAsPercent",

  "nIntersectSignifGOid_manyAsTopTADs" = "# signif. GO (id) - intersect manyAsTopTADs",
  "nIntersectSignifGOterm_manyAsTopTADs" = "# signif. GO (term) - intersect manyAsTopTADs",
  "nIntersectSignifGOmin_manyAsTopTADs" = "# signif. GO (min. set) - intersect manyAsTopTADs",
  "intersectSignifGOidRatio_manyAsTopTADs" = "ratio signif. GO (id) - intersect manyAsTopTADs",
  "intersectSignifGOtermRatio_manyAsTopTADs" = "ratio signif. GO (term) - intersect manyAsTopTADs",
  "intersectSignifGOminRatio_manyAsTopTADs" = "ratio signif. GO (min. set) - intersect manyAsTopTADs",


  "aucFCC" = "AUC ratio - FCC",
  "aucCoexprDist" = "AUC ratio - coexpr.",

  "nTopTADs_genes_signifGOidRatio" = "ratio # signif. GO (id) - topTADs_genes",
  "nTopTADs_genes_signifGOtermRatio" = "ratio # signif. GO (term) - topTADs_genes",
  "nTopTADs_genes_signifGOminRatio" = "ratio # signif. GO (min. set) - topTADs_genes",

  "nTopGenes_manyAsPercent_signifGOidRatio" = "ratio # signif. GO (id) - topGenes_manyAsPercent",
  "nTopGenes_manyAsPercent_signifGOtermRatio" = "ratio # signif. GO (term) - topGenes_manyAsPercent",
  "nTopGenes_manyAsPercent_signifGOminRatio" = "ratio # signif. GO (min. set) - topGenes_manyAsPercent",

  "nTopGenes_manyAsTopTADs_signifGOidRatio" = "ratio # signif. GO (id) - topGenes_manyAsTopTADs",
  "nTopGenes_manyAsTopTADs_signifGOtermRatio" = "ratio # signif. GO (term) - topGenes_manyAsTopTADs",
  "nTopGenes_manyAsTopTADs_signifGOminRatio" = "ratio # signif. GO (min. set) - topGenes_manyAsTopTADs",

  "topGenes_manyAsPercent_intersectRatio" = "ratio # signif. genes in intersect - topGenes_manyAsPercent",

  "topGenes_manyAsTopTADs_intersectRatio" = "ratio # signif. genes in intersect - topGenes_manyAsTopTADs",

  "topTADs_genes_intersectRatio_manyAsPercent" = "ratio # signif. genes in intersect - topTADs_genes",
  "nIntersectGenes_manyAsPercent" = "# signif. genes - intersect manyAsPercent",
  "nUnionGenes_manyAsPercent" = "# signif. genes - union",
  "intersectGenesRatio_manyAsPercent" = "ratio signif. genes - intersect manyAsPercent",

  "topTADs_genes_intersectRatio_manyAsTopTADs" = "ratio # signif. genes in intersect - topTADs_genes",
  "nIntersectGenes_manyAsTopTADs" = "# signif. genes - intersect manyAsTopTADs",
  "nUnionGenes_manyAsTopTADs" = "# signif. genes - union",
  "intersectGenesRatio_manyAsTopTADs" = "ratio signif. genes - intersect manyAsTopTADs",

"topGenes_manyAsTopTADs_meanGOdim" = "Avg. ratio # GO genes from list/# GO genes - topGenes_manyAsTopTADs",
"topGenes_manyAsPercent_meanGOdim" = "Avg. ratio # GO genes from list/# GO genes - topGenes_manyAsPercent",
"topTADs_genes_meanGOdim" = "Avg. ratio # GO genes from list/# GO genes - topTADs_genes"
)

GO_aliases_top <- gsub(" - topGenes_manyAsPercent", "", GO_aliases_top)
GO_aliases_top <- gsub(" - topTADs_genes", "", GO_aliases_top)


GO_offSets_top <- c(
  "dataset" = NA,
  "nTADs" = 500,
  "nTopTADs" = 100, 
  "nTopTADs_genes" = 100,
  "nTopTADs_genes_signifGOid" = 50,
  "nTopTADs_genes_signifGOterm" = 50,
  "nTopTADs_genes_signifGOmin" = 50,
  "topTADs_genes_g_density" = 0.01,
  "topTADs_genes_g_diameter" = 1,
  "topTADs_genes_g_meanDist" = 0.1,
  "topTADs_genes_g_meanDistDir" = 0.1,
  "topTADs_genes_g_meanEccentricity" = 0.1,
  "topTADs_genes_g_meanBetweenness" = 1,

  "nGenes" = 1000,
  "nTopGenes_manyAsPercent" = 100, 
  "nTopGenes_manyAsPercent_signifGOid" = 50,
  "nTopGenes_manyAsPercent_signifGOterm" = 50,
  "nTopGenes_manyAsPercent_signifGOmin" = 50,
  "topGenes_manyAsPercent_g_density" = 0.01,
  "topGenes_manyAsPercent_g_diameter"= 1,
  "topGenes_manyAsPercent_g_meanDist" = 0.1, 
  "topGenes_manyAsPercent_g_meanDistDir" = 0.1, 
  "topGenes_manyAsPercent_g_meanEccentricity" = 0.1, 
  "topGenes_manyAsPercent_g_meanBetweenness" = 1,

  "nTopGenes_manyAsTopTADs" = 100, 
  "nTopGenes_manyAsTopTADs_signifGOid" = 50,
  "nTopGenes_manyAsTopTADs_signifGOterm" = 50,
  "nTopGenes_manyAsTopTADs_signifGOmin" = 50,
  "topGenes_manyAsTopTADs_g_density" = 0.01,
  "topGenes_manyAsTopTADs_g_diameter"= 1,
  "topGenes_manyAsTopTADs_g_meanDist" = 0.1, 
  "topGenes_manyAsTopTADs_g_meanDistDir" = 0.1, 
  "topGenes_manyAsTopTADs_g_meanEccentricity" = 0.1, 
  "topGenes_manyAsTopTADs_g_meanBetweenness" = 1,

  "nIntersectSignifGOid_manyAsPercent" = 1,
  "nIntersectSignifGOterm_manyAsPercent" = 1,
  "nIntersectSignifGOmin_manyAsPercent" = 1,
  "intersectSignifGOidRatio_manyAsPercent" = 0.01,
  "intersectSignifGOtermRatio_manyAsPercent" = 0.01,
  "intersectSignifGOminRatio_manyAsPercent" = 0.01,

  "nIntersectSignifGOid_manyAsTopTADs" = 1,
  "nIntersectSignifGOterm_manyAsTopTADs" = 1,
  "nIntersectSignifGOmin_manyAsTopTADs" = 1,
  "intersectSignifGOidRatio_manyAsTopTADs" = 0.01,
  "intersectSignifGOtermRatio_manyAsTopTADs" = 0.01,
  "intersectSignifGOminRatio_manyAsTopTADs" = 0.01,

  "aucFCC" = 0.03,
  "aucCoexprDist" = 0.03,

  "nTopTADs_genes_signifGOidRatio" = 0.001,
  "nTopTADs_genes_signifGOtermRatio" = 0.001,
  "nTopTADs_genes_signifGOminRatio" = 0.001,

  "nTopGenes_manyAsPercent_signifGOidRatio" = 0.001,
  "nTopGenes_manyAsPercent_signifGOtermRatio" = 0.001,
  "nTopGenes_manyAsPercent_signifGOminRatio" = 0.001,

  "nTopGenes_manyAsTopTADs_signifGOidRatio" = 0.001,
  "nTopGenes_manyAsTopTADs_signifGOtermRatio" = 0.001,
  "nTopGenes_manyAsTopTADs_signifGOminRatio" = 0.001,
    
  "topGenes_manyAsTopTADs_intersectRatio" = 0.01,    
  "topGenes_manyAsPercent_intersectRatio" = 0.01,

  "topTADs_genes_intersectRatio_manyAsTopTADs" = 0.01,
  "topTADs_genes_intersectRatio_manyAsPercent" = 0.01,

  "nIntersectGenes_manyAsPercent" = 1,
  "nUnionGenes_manyAsPercent" = 1,
  "nIntersectGenes_manyAsTopTADs" = 1,
  "nUnionGenes_manyAsTopTADs" = 1,

  "intersectGenesRatio_manyAsTopTADs" = 0.01,
  "intersectGenesRatio_manyAsPercent" = 0.01,


"topGenes_manyAsTopTADs_meanGOdim" = 0.01,
"topGenes_manyAsPercent_meanGOdim" = 0.01,
"topTADs_genes_meanGOdim" = 0.01

)


GO_aliases_common_top <- c(
    "signifGOid"="# signif. GO (id)",
    "signifGOterm"="# signif. GO (term)",
    "signifGOmin"="# signif. GO (min. set)",
    "g_density"="graph density",
    "g_diameter"="graph diameter",
    "g_meanDist"="graph mean dist. (undir.)",
    "g_meanDistDir"="graph mean dist. (dir.)",
    "g_meanEccentricity"="graph mean eccentricity",
    "g_meanBetweenness"="graph mean betweenness",
    "signifGOidRatio"="ratio # signif. GO (id)",
    "signifGOtermRatio"="ratio # signif. GO (term)",
    "signifGOminRatio"="ratio # signif. GO (min. set)",
    "nTop"="# selected features",
    "intersectRatio" = "ratio # intersect genes",
    "meanGOdim" = "Avg. ratio # GO genes from list/# GO genes"



)                         


GO_legPos_top <- rep("bottomright", length(GO_aliases_top) )
names(GO_legPos_top) <- names(GO_aliases_top)

topRight_var <- c(
"topGenes_manyAsPercent_intersectRatio",
"topGenes_manyAsTopTADs_intersectRatio",
"topTADs_genes_intersectRatio",
"intersectGenesRatio",

"intersectSignifGOtermRatio",
"nTopGenes_manyAsPercent_signifGOtermRatio",
"nTopGenes_manyAsTopTADs_signifGOtermRatio",
"nTopTADs_genes_signifGOtermRatio"
)
botLeft_var <- c(
"topTADs_genes_g_meanEccentricity",
"topGenes_manyAsPercent_g_meanEccentricity",
"topGenes_manyAsTopTADs_g_meanEccentricity",

"topTADs_genes_g_meanDist",
"topGenes_manyAsPercent_g_meanDist",
"topGenes_manyAsTopTADs_g_meanDist",

"topTADs_genes_g_meanDistDir",
"topGenes_manyAsPercent_g_meanDistDir",
"topGenes_manyAsTopTADs_g_meanDistDir",

"topTADs_genes_g_diameter",
"topGenes_manyAsPercent_g_diameter",
"topGenes_manyAsTopTADs_g_diameter",

"topTADs_genes_g_density",
"topGenes_manyAsPercent_g_density",
"topGenes_manyAsTopTADs_g_density",

"topTADs_genes_g_meanBetweenness",
"topGenes_manyAsPercent_g_meanBetweenness",
"topGenes_manyAsTopTADs_g_meanBetweenness",

"intersectSignifGOminRatio",
"nTopGenes_manyAsPercent_signifGOminRatio",
"nTopGenes_manyAsTopTADs_signifGOminRatio",
"nTopTADs_genes_signifGOminRatio"
)


GO_legPos_top[names(GO_legPos_top) %in% topRight_var ] <- "topright"
GO_legPos_top[names(GO_legPos_top) %in% botLeft_var ] <- "bottomleft"

######################################################################################################################################################################################################
###################################################################################################################################################################################################### MAPPING ALIASES FOR GO VARIABLES (pvalSelect) (vector)
######################################################################################################################################################################################################
GO_aliases_pvalSelect <- c(
  "dataset" = "dataset",
  "nTADs" = "# TADs",
  "nSelectTADs" = "# selectTADs",
  "nSelectTADs_genes" = "# selectTADs_genes",        
  "nSelectTADs_genes_signifGOid" = "# signif. GO (id) - selectTADs_genes",
  "nSelectTADs_genes_signifGOterm" = "# signif. GO (term) - selectTADs_genes",
  "nSelectTADs_genes_signifGOmin" = "# signif. GO (min. set) - selectTADs_genes",
  "selectTADs_genes_g_density" = "graph density - selectTADs_genes",
  "selectTADs_genes_g_diameter" = "graph diameter - selectTADs_genes",
  "selectTADs_genes_g_meanDist" = "graph mean dist. (undir.) - selectTADs_genes",
  "selectTADs_genes_g_meanDistDir" = "graph mean dist. (dir.) - selectTADs_genes",
  "selectTADs_genes_g_meanEccentricity" = "graph mean eccentricity - selectTADs_genes",
  "selectTADs_genes_g_meanBetweenness" = "graph mean betweenness - selectTADs_genes",
  
  "nGenes" = "# genes",
  "nSelectGenes" = "# selectGenes",
  "nSelectGenes_signifGOid" = "# signif. GO (id) - selectGenes",
  "nSelectGenes_signifGOterm" = "# signif. GO (term) - selectGenes",
  "nSelectGenes_signifGOmin" = "# signif. GO (min. set) - selectGenes",
  "selectGenes_g_density" = "graph density - selectGenes",
  "selectGenes_g_diameter"= "graph diameter - selectGenes",
  "selectGenes_g_meanDist" = "graph mean dist. (undir.) - selectGenes",
  "selectGenes_g_meanDistDir" = "graph mean dist. (dir.) - selectGenes",
  "selectGenes_g_meanEccentricity" = "graph mean eccentricity - selectGenes",
  "selectGenes_g_meanBetweenness" ="graph mean betweenness - selectGenes",
  "nIntersectSignifGOid" = "# signif. GO (id) - intersect",
  "nIntersectSignifGOterm" = "# signif. GO (term) - intersect",
  "nIntersectSignifGOmin" = "# signif. GO (min. set) - intersect",
  "intersectSignifGOidRatio" = "ratio signif. GO (id) - intersect",
  "intersectSignifGOtermRatio" = "ratio signif. GO (term) - intersect",
  "intersectSignifGOminRatio" = "ratio signif. GO (min. set) - intersect",
  "aucFCC" = "AUC ratio - FCC",
  "aucCoexprDist" = "AUC ratio - coexpr.",
  "nSelectTADs_genes_signifGOidRatio" = "ratio # signif. GO (id) - selectTADs_genes",
  "nSelectTADs_genes_signifGOtermRatio" = "ratio # signif. GO (term) - selectTADs_genes",
  "nSelectTADs_genes_signifGOminRatio" = "ratio # signif. GO (min. set) - selectTADs_genes",
  "nSelectGenes_signifGOidRatio" = "ratio # signif. GO (id) - selectGenes",
  "nSelectGenes_signifGOtermRatio" = "ratio # signif. GO (term) - selectGenes",
  "nSelectGenes_signifGOminRatio" = "ratio # signif. GO (min. set) - selectGenes",

  "selectGenes_intersectRatio" = "ratio # signif. genes in intersect - selectGenes",
  "selectTADs_genes_intersectRatio" = "ratio # signif. genes in intersect - selectTADs_genes",
  "nIntersectGenes" = "# signif. genes - intersect",
  "nUnionGenes" = "# signif. genes - union",
  "intersectGenesRatio" = "ratio signif. genes - intersect",


"selectGenes_meanGOdim" = "Avg. ratio # GO genes from list/# GO genes - selectGenes",
"selectTADs_genes_meanGOdim" = "Avg. ratio # GO genes from list/# GO genes - selectTADs_genes"

)

GO_aliases_pvalSelect <- gsub(" - selectGenes", "", GO_aliases_pvalSelect)
GO_aliases_pvalSelect <- gsub(" - selectTADs_genes", "", GO_aliases_pvalSelect)


GO_offSets_pvalSelect <- c(
  "dataset" = NA,
  "nTADs" = 500,
  "nSelectTADs" = 100, 
  "nSelectTADs_genes" = 100,
  "nSelectTADs_genes_signifGOid" = 50,
  "nSelectTADs_genes_signifGOterm" = 50,
  "nSelectTADs_genes_signifGOmin" = 50,
  "selectTADs_genes_g_density" = 0.01,
  "selectTADs_genes_g_diameter" = 1,
  "selectTADs_genes_g_meanDist" = 0.1,
  "selectTADs_genes_g_meanDistDir" = 0.1,
  "selectTADs_genes_g_meanEccentricity" = 0.1,
  "selectTADs_genes_g_meanBetweenness" = 1,
  "nGenes" = 1000,
  "nSelectGenes" = 100, 
  "nSelectGenes_signifGOid" = 50,
  "nSelectGenes_signifGOterm" = 50,
  "nSelectGenes_signifGOmin" = 50,
  "selectGenes_g_density" = 0.01,
  "selectGenes_g_diameter"= 1,
  "selectGenes_g_meanDist" = 0.1, 
  "selectGenes_g_meanDistDir" = 0.1, 
  "selectGenes_g_meanEccentricity" = 0.1, 
  "selectGenes_g_meanBetweenness" = 1,
  "nIntersectSignifGOid" = 1,
  "nIntersectSignifGOterm" = 1,
  "nIntersectSignifGOmin" = 1,
  "intersectSignifGOidRatio" = 0.01,
  "intersectSignifGOtermRatio" = 0.01,
  "intersectSignifGOminRatio" = 0.01,
  "aucFCC" = 0.03,
  "aucCoexprDist" = 0.03,
  "nSelectTADs_genes_signifGOidRatio" = 0.001,
  "nSelectTADs_genes_signifGOtermRatio" = 0.001,
  "nSelectTADs_genes_signifGOminRatio" = 0.001,
  "nSelectGenes_signifGOidRatio" = 0.001,
  "nSelectGenes_signifGOtermRatio" = 0.001,
  "nSelectGenes_signifGOminRatio" = 0.001,
    
  "selectGenes_intersectRatio" = 0.01,
  "selectTADs_genes_intersectRatio" = 0.01,
  "nIntersectGenes" = 1,
  "nUnionGenes" = 1,
  "intersectGenesRatio" = 0.01,

"selectGenes_meanGOdim" = 0.01,
"selectTADs_genes_meanGOdim" = 0.01



)


GO_aliases_common_pvalSelect <- c(
    "signifGOid"="# signif. GO (id)",
    "signifGOterm"="# signif. GO (term)",
    "signifGOmin"="# signif. GO (min. set)",
    "g_density"="graph density",
    "g_diameter"="graph diameter",
    "g_meanDist"="graph mean dist. (undir.)",
    "g_meanDistDir"="graph mean dist. (dir.)",
    "g_meanEccentricity"="graph mean eccentricity",
    "g_meanBetweenness"="graph mean betweenness",
    "signifGOidRatio"="ratio # signif. GO (id)",
    "signifGOtermRatio"="ratio # signif. GO (term)",
    "signifGOminRatio"="ratio # signif. GO (min. set)",
    "nSelect"="# selected features",
    "intersectRatio" = "ratio # intersect genes",
"meanGOdim" = "Avg. ratio # GO genes from list/# GO genes"
)                         


GO_legPos_pvalSelect <- rep("bottomright", length(GO_aliases_pvalSelect) )
names(GO_legPos_pvalSelect) <- names(GO_aliases_pvalSelect)

topRight_var <- c(
"selectGenes_intersectRatio",
"selectTADs_genes_intersectRatio",
"intersectGenesRatio",

"intersectSignifGOtermRatio",
"nSelectGenes_signifGOtermRatio",
"nSelectTADs_genes_signifGOtermRatio"
)
botLeft_var <- c(
"selectTADs_genes_g_meanEccentricity",
"selectGenes_g_meanEccentricity",

"selectTADs_genes_g_meanDist",
"selectGenes_g_meanDist",

"selectTADs_genes_g_meanDistDir",
"selectGenes_g_meanDistDir",

"selectTADs_genes_g_diameter",
"selectGenes_g_diameter",

"selectTADs_genes_g_density",
"selectGenes_g_density",

"selectTADs_genes_g_meanBetweenness",
"selectGenes_g_meanBetweenness",

"intersectSignifGOminRatio",
"nSelectGenes_signifGOminRatio",
"nSelectTADs_genes_signifGOminRatio"
)



GO_legPos_pvalSelect[names(GO_legPos_pvalSelect) %in% topRight_var ] <- "topright"
GO_legPos_pvalSelect[names(GO_legPos_pvalSelect) %in% botLeft_var ] <- "bottomleft"





######################################################################################################################################################################################################
###################################################################################################################################################################################################### draw_myVenn (function)
######################################################################################################################################################################################################


draw_myVenn <- function(var1, var2, varIntersect, tit, subTit, dataDT, subPatt1="", subPatt2="") {
  
  stopifnot(var1 %in% colnames(dataDT))
  stopifnot(var2 %in% colnames(dataDT))
  stopifnot(varIntersect %in% colnames(dataDT))
  
  var1_cat <- gsub(subPatt1, "", var1)
  var2_cat <- gsub(subPatt2, "", var2)
  
  a1 <- dataDT[1, var1]
  a2 <- dataDT[1, var2]
  aIntersect <- dataDT[1, varIntersect]
  
  a1 <- ifelse(is.na(a1), 0, a1)
  a2 <- ifelse(is.na(a2), 0, a2)
  aIntersect <- ifelse(is.na(aIntersect), 0, aIntersect)
  
#  cat("a1 = ", a1, "\n")
#  cat("a2 = ", a2, "\n")
#  cat("aIntersect = ", aIntersect, "\n")
  
  
  grid.newpage()

  if(all(c(a1, a2, aIntersect) == 0)){
  vD <- textGrob(paste0("(! all 0 !)"),gp=gpar(fontsize=16, fontface="italic"))
  } else {
  vD <- draw.pairwise.venn(area1=a1, area2=a2, cross.area=aIntersect, 
                           category = c(var1_cat, var2_cat), # => this will draw 11-11-9
                           lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                           cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
  }
  # titGrob <- textGrob(paste0(tit),gp=gpar(fontsize=20, fontface="bold"))
  subGrob <- textGrob(paste0(subTit),gp=gpar(fontsize=16, fontface="bold"))
  # retVD <- grid.arrange(gTree(children=vD), top=titGrob, bottom=subGrob)

  if(all(c(a1, a2, aIntersect) == 0)){
    retVD <- grid.arrange(vD, bottom=subGrob)
  } else {
   retVD <- grid.arrange(gTree(children=vD), bottom=subGrob)
  }
  return(retVD)
}


######################################################################################################################################################################################################
###################################################################################################################################################################################################### not_used_yet_plot_cumsumDiff05_withLines (function)
######################################################################################################################################################################################################

not_used_yet_plot_cumsumDiff05_withLines <- function(observ_vect, permut_DT, 
           all_quantThresh = NULL,
		   pointObsCol = "black", 
           my_stat = "ratioDown",
           my_main = NULL,
           my_ylab = NULL,
           my_xlab = NULL,
           cexMain = 1,
           polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),			
		   linePermutCol = "navyblue",
           departureValue = 0.5, drawline=FALSE) {
  
  plotTrueType <- ifelse(drawline, "l", "p")

  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)

  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))


  observed_yval <- cumsum(abs(observ_vect - departureValue))

  permut_yval_min <- apply(diff_05_permut, 1, function(x) min(x))
  permut_yval_max <- apply(diff_05_permut, 1, function(x) max(x))
  permut_yval_mean <- apply(diff_05_permut, 1, function(x) mean(x))
  stopifnot(!any(is.na(permut_yval_min)))
  stopifnot(!any(is.na(permut_yval_max)))
  stopifnot(!any(is.na(permut_yval_mean)))
  
  plot(observed_yval ~ x_val,
       main= my_main,
       cex.main = cexMain,
       type = plotTrueType,
       pch = 16, cex = 0.7,
       xlab= my_xlab, 
       ylab= my_ylab,
       bty="l")
  polygon(x = c(x_val, rev(x_val)), 
          y = c( permut_yval_min, rev(permut_yval_max)),
          border=NA,
          col = polygonPermutCol)

  # add lines:
  lines(x = c(x_val),
      y = c(permut_yval_min),
      col = linePermutCol)

  lines(x = c(x_val),
      y = c(permut_yval_max),
      col = linePermutCol)

  lines(x = c(x_val),
      y = c(permut_yval_mean),
      col = linePermutCol)

  if(!is.null(all_quantThresh)) {
	stopifnot(is.numeric(all_quantThresh))
    stopifnot(all_quantThresh >= 0 & all_quantThresh <= 1)

    all_auc_quantThresh <- list()

	for(quantThresh in all_quantThresh) {
		permut_yval_currThresh <- apply(diff_05_permut, 1, function(x) quantile(x, probs = quantThresh, names=F) )

	    lines(x = c(x_val),
		      y = c(permut_yval_currThresh),

		      col = linePermutCol)

         aucPermut_currThresh <- auc(x = x_val, y = permut_yval_currThresh)
		 all_auc_quantThresh[[as.character(quantThresh)]] <- aucPermut_currThresh

	}

    legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15, -1),
         lty = c(-1, -1, 1), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut", paste0("min-mean-max-quant. thresh: ", paste0(all_quantThresh, collapse=", "))), 
         pt.cex = c(0.7, 2, -1),
         col = c(pointObsCol, polygonPermutCol, linePermutCol),
         bty="n")


    all_auc_quantThresh_vect <- setNames(unlist(all_auc_quantThresh), paste0("quant_", names(all_auc_quantThresh)))

  } else {

    legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15, -1),
         lty = c(-1, -1, 1), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut", paste0("min-mean-max")), 
         pt.cex = c(0.7, 2, -1),
         col = c(pointObsCol, polygonPermutCol, linePermutCol),
         bty="n")

	all_auc_quantThresh_vect <- NULL

}


  # compute and return the AUC

  auc_values <- c(
	observed_auc = auc(x = x_val, y = observed_yval),
	minPermut_auc = auc(x = x_val, y = permut_yval_min),
	meanPermut_auc = auc(x = x_val, y = permut_yval_mean),
	maxPermut_auc = auc(x = x_val, y = permut_yval_max)
	)
  auc_values <- c(auc_values, all_auc_quantThresh_vect)

	return(auc_values)

}

