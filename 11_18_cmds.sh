#!/usr/bin/bash
# 
# ./11_18_cmds.sh

start_time=$(date -R)    


#Rscript aucFCC_vs_aucCoexprDist_cancerDS.R

#Rscript aucFCC_vs_aucCoexprDist_TCGA.R

Rscript aucFCC_vs_aucCoexprDist.R

############################################

#Rscript gene_variance.R log2fpkm
#Rscript gene_variance_cancerDS.R log2fpkm
#Rscript gene_variance_TCGA.R log2fpkm
#Rscript gene_variance_noTCGA.R log2fpkm

############################################

# Rscript cumul_gene_variance.R log2fpkm 1000 raw
# Rscript cumul_gene_variance.R log2fpkm 1000 rescVar
# Rscript cumul_gene_variance.R log2fpkm 1000 log10

# Rscript cumul_gene_variance.R log2fpkm 1000 raw_cropVar
# Rscript cumul_gene_variance.R log2fpkm 1000 rescVar_cropVar
# Rscript cumul_gene_variance.R log2fpkm 1000 log10_cropVar


# Rscript cumul_gene_FC.R 1000 raw_cropFC
# Rscript cumul_gene_FC.R 1000 rescFC_cropFC
# Rscript cumul_gene_FC.R 1000 log10_cropFC

# Rscript cumul_gene_FC.R 1000 raw
# Rscript cumul_gene_FC.R 1000 rescFC
# Rscript cumul_gene_FC.R 1000 log10

############################################

#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 5 avg
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 10 avg
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 15 avg
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 20 avg
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 25 avg
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 66 avg

#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 5 FCC
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 10 FCC
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 15 FCC
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 20 FCC
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 25 FCC
#Rscript cmp_DE_TADs_genes_GO_pvalSelect_top.R 66 FCC

###############################################

#Rscript cmp_DE_TADs_genes_GO_top.R 5 avg
#Rscript cmp_DE_TADs_genes_GO_top.R 10 avg
#Rscript cmp_DE_TADs_genes_GO_top.R 15 avg
#Rscript cmp_DE_TADs_genes_GO_top.R 20 avg
#Rscript cmp_DE_TADs_genes_GO_top.R 25 avg
#Rscript cmp_DE_TADs_genes_GO_top.R 66 avg

#Rscript cmp_DE_TADs_genes_GO_top.R 5 FCC
#Rscript cmp_DE_TADs_genes_GO_top.R 10 FCC
#Rscript cmp_DE_TADs_genes_GO_top.R 15 FCC
#Rscript cmp_DE_TADs_genes_GO_top.R 20 FCC
#Rscript cmp_DE_TADs_genes_GO_top.R 25 FCC
#Rscript cmp_DE_TADs_genes_GO_top.R 66 FCC


#Rscript AUC_ratio_cumulVect.R FC log10
#Rscript AUC_ratio_cumulVect.R FC log10_cropFC
#Rscript AUC_ratio_cumulVect.R FC raw
#Rscript AUC_ratio_cumulVect.R FC raw_cropFC
#Rscript AUC_ratio_cumulVect.R FC rescFC
#Rscript AUC_ratio_cumulVect.R FC rescFC_cropFC

#Rscript AUC_ratio_cumulVect.R var log10
#Rscript AUC_ratio_cumulVect.R var log10_cropVar
#Rscript AUC_ratio_cumulVect.R var raw
#Rscript AUC_ratio_cumulVect.R var raw_cropVar
#Rscript AUC_ratio_cumulVect.R var rescVar
#Rscript AUC_ratio_cumulVect.R var rescVar_cropVar



###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0



