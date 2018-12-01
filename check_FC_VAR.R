
setDir <- "/media/electron"
load("")

load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18", "CUMUL_GENE_FC/1000/raw/all_ds_geneFC_DT.Rdata"))
all_ds_geneFC_DT <- all_ds_geneFC_DT[order(all_ds_geneFC_DT$absFC, decreasing = T),]
head(all_ds_geneFC_DT)


unique(all_ds_geneFC_DT$dataset[all_ds_geneFC_DT$absFC > 1000])
# "GSE77509_normal_tumor"    "GSE77509_normal_ptt"      "GSE77509_ptt_tumor"       "GSE68719_norm_park"       "GSE64810_control_carrier"
# [6] "GSE101521_control_mdd"

unique(all_ds_geneFC_DT$dataset[all_ds_geneFC_DT$absFC > 50])
# [1] "GSE77509_normal_tumor"    "GSE77509_normal_ptt"      "GSE77509_ptt_tumor"       "GSE68719_norm_park"       "GSE64810_control_carrier"
# [6] "GSE101521_control_mdd"    "GSE77314_normal_tumor"   

load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18", "CUMUL_GENE_VARIANCE/LOG2FPKM_1000/raw/all_ds_geneVarDT.Rdata"))
all_ds_geneVarDT <- all_ds_geneVarDT[order(all_ds_geneVarDT$var, decreasing = T),]
head(all_ds_geneVarDT)
unique(all_ds_geneVarDT$dataset[all_ds_geneVarDT$var > 1000])



x=load(paste0(setDir, "/mnt/ed4/marie/other_datasets/GSE77509/20.09_prepData/GSE77509_count_all_DT.Rdata"))
s1=load(paste0(setDir, "/mnt/ed4/marie/other_datasets/GSE77509/20.09_prepData/normal_ID.Rdata"))                                                                                                                                                                                     
s2=load(paste0(setDir, "/mnt/ed4/marie/other_datasets/GSE77509/20.09_prepData/tumor_ID.Rdata"))                                                                                                                                                                                                
s1
s2
# [1] "tumor_ID"
x
# [1] "GSE77509_count_all_DT"
gene ="ENSG00000163631"                                                                                                                                                                                                                                                          
dt1=GSE77509_count_all_DT[gene,s1]
dt1
NULL
dt1=GSE77509_count_all_DT[gene,s1]
s1
# [1] "normal_ID"
dt1=GSE77509_count_all_DT[gene,normal_ID]                                                                                                                                                                                                                                        
dt1
#                      N3      N6       N7      N8     N10     N11     N12
# ENSG00000163631 2736739 6047006 13257400 7626439 4342366 5960678 4425673
#                     N13     N14      N15     N16      N17     N18      N19
# ENSG00000163631 2893864 8568325 16092854 7405320 11651223 9235849 14807387
#                      N20     N21     N22      N24      N25      N26
# ENSG00000163631 10956388 9896547 8451162 12359975 10612323 10298930
dt2=GSE77509_count_all_DT[gene,tumor_ID]                                                                                                                                                                                                                                         
dt2
#                      T3      T6      T7      T8      T10     T11     T12
# ENSG00000163631 2746327 2373368 1498646 1747692 12742.96 1122509 2283326
#                      T13      T14     T15      T16     T17     T18     T19
# ENSG00000163631 457071.4 344016.5 3376222 38413.22 1453831 5453040 1301789
#                     T20     T21     T22     T24    T25      T26
# ENSG00000163631 3915292 3780290 1340872 2971422 642548 158262.2
rowMeans(dt1)
# ENSG00000163631 
#         8881322 
rowMeans(dt2)                                                                                                                                                                                                                                                                    
# ENSG00000163631 
#         1850884 
rowMeans(dt1)/rowMeans(dt2)
# ENSG00000163631 
#        4.798422 

x=load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/GSE77509_normal_tumor/0_prepGeneData/rna_rnaseqDT.Rdata"))
x
# [1] "rna_rnaseqDT"
rt1=rna_rnaseqDT[gene,normal_ID]                                                                                                                                                                                                                                                 
rt1
#                      N3      N6       N7      N8     N10     N11     N12
# ENSG00000163631 2736739 6047006 13257400 7626439 4342366 5960678 4425673
#                     N13     N14      N15     N16      N17     N18      N19
# ENSG00000163631 2893864 8568325 16092854 7405320 11651223 9235849 14807387
#                      N20     N21     N22      N24      N25      N26
# ENSG00000163631 10956388 9896547 8451162 12359975 10612323 10298930
rt2=rna_rnaseqDT[gene,tumor_ID]                                                                                                                                                                                                                                                  
rt2
#                      T3      T6      T7      T8      T10     T11     T12
# ENSG00000163631 2746327 2373368 1498646 1747692 12742.96 1122509 2283326
#                      T13      T14     T15      T16     T17     T18     T19
# ENSG00000163631 457071.4 344016.5 3376222 38413.22 1453831 5453040 1301789
#                     T20     T21     T22     T24    T25      T26
# ENSG00000163631 3915292 3780290 1340872 2971422 642548 158262.2
rowMeans(rt1)/rowMeans(rt2)
# ENSG00000163631 
       # 4.798422 

