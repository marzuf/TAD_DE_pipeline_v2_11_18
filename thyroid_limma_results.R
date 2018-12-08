setDir <- "/media/electron"
deTableFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER/TCGAthca_mutBRAF_mutKRAS", "1_runGeneDE", "DE_topTable.Rdata")

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gffDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
gffDT$entrezID <- as.character(gffDT$entrezID)
  
  
DE_topTableDT <- eval(parse(text=load(deTableFile)))
DE_topTableDT$genes <- as.character(DE_topTableDT$genes)
DE_topTableDT <- DE_topTableDT[DE_topTableDT$genes %in% gffDT$entrezID,]
DE_topTableDT$symbols <- unlist(sapply(as.character(DE_topTable$genes), function(x) gffDT$symbol[gffDT$entrezID == x]))
DE_topTableDT <- DE_topTableDT[order(abs(DE_topTableDT$logFC), decreasing = T),]
