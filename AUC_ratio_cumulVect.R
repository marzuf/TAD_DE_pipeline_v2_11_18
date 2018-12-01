library(foreach)
library(doMC)
library(flux)
library(ggplot2)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

startTime <- Sys.time()

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

barCol <- "blue4"

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)

source("analysis_utils.R")

# Rscript AUC_ratio_cumulVect.R <cumulVar>
# Rscript AUC_ratio_cumulVect.R FC log10

nTop <- 1000 # how many top features were used to build the data

cat("> START: AUC_ratio_cumulVect.R\n")

cumulVar <- "FC"
dataType <- "log10"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
cumulVar <- args[1]
dataType <- args[2]
stopifnot(cumulVar %in% c("FC", "var"))
stopifnot(dataType %in% c("log10", "raw", "rescVar", "rescFC", "log10_cropFC", "log10_cropVar", "raw_cropFC", "raw_cropVar", "rescVar_cropVar", "rescFC", "rescFC_cropFC"))

outFold <- file.path(paste0("AUC_RATIO_CUMULVECT"), cumulVar, dataType)
system(paste0("mkdir -p ", outFold))

if(cumulVar == "FC") {
  
  vectFile <- file.path("CUMUL_GENE_FC",
                        nTop,
                        dataType,
                        paste0("all_y_", dataType, ".Rdata"))
  
  cat("vectFile = ", vectFile, "\n")
  
  stopifnot(file.exists(vectFile))
  
  all_yvals <- eval(parse(text = load(vectFile)))
  all_yvals <- all_yvals[["mostFC"]]
  cumulVarTit <- "Gene expression FC"
} else if(cumulVar == "var") {
  
  vectFile <- file.path("CUMUL_GENE_VARIANCE",
                        paste0("LOG2FPKM_", nTop),
                        dataType,
                        paste0("all_y_", dataType, ".Rdata"))
  
  stopifnot(file.exists(vectFile))
  
  all_yvals <- eval(parse(text = load(vectFile)))
  all_yvals <- all_yvals[["mostVar"]]
  cumulVarTit <- "Gene expression variance" 
}


# COMPUTE THE AUC

auc_values <- lapply(all_yvals, function(x) {
  auc(x = c(1:length(x)),
            y = x)  
})

maxAUC <- max(unlist(auc_values))

resc_auc_values <- lapply(auc_values, function(x) x/maxAUC)

stopifnot(unlist(resc_auc_values) >= 0 & unlist(resc_auc_values) <= 1)


auc_values <- list(AUC = auc_values,
                   rescAUC = resc_auc_values)
outFile <- file.path(outFold, "auc_values.Rdata")
save(auc_values, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################
########################################################################################## RETRIEVE FCC AUC RATIO 
##########################################################################################

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)

all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  if(!file.exists(aucFCC_file)) cat("... aucFCC_file = ", aucFCC_file, "\n")
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

all_ds_aucFCC <- sort(all_ds_aucFCC, decreasing=TRUE)

stopifnot(names(resc_auc_values) %in% names(all_ds_aucFCC))


######################################################################################
###################################################################################### AUC values vs FCC
######################################################################################

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

cat("!!! TEMPORARY FILTER DS !!!")
all_ds <- all_ds[all_ds %in% names(resc_auc_values)]

resc_auc <- setNames(unlist(resc_auc_values),names(resc_auc_values))
stopifnot(names(resc_auc) %in% all_ds)
stopifnot(all_ds %in% names(resc_auc))

resc_auc <- resc_auc[all_ds]
all_ds_aucFCC <- all_ds_aucFCC[all_ds]
stopifnot(names(resc_auc) == names(all_ds_aucFCC))

subTit <- paste0("all datasets (n=", length(all_ds), ")")

stopifnot(all_ds %in% names(dataset_proc_colors))
curr_colors <- dataset_proc_colors[all_ds]

outFile <- file.path(outFold, paste0(cumulVar, "_vs_FCC.", plotType))

all_ds_aucFCC <- (all_ds_aucFCC-1)*100
resc_auc <- (resc_auc)*100

offsetRange <- c(-15,15)

#myTit <- paste0("% increase AUC ", cumulVar, " vs. FCC")
#myxlab <- paste0("% increase AUC - FCC")
#myylab <- paste0("% increase AUC - ", cumulVar, " (", dataType, ")" )

myTit <- paste0("% rel. AUC ", cumulVar, " vs. % increase AUC FCC")
myxlab <- paste0("% increase AUC - FCC")
myylab <- paste0("% rel. AUC - ", cumulVar, " (", dataType, ")" )


do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(x = all_ds_aucFCC,
     y = resc_auc,
     xlim = range(all_ds_aucFCC) + offsetRange,
     ylim = range(resc_auc) + offsetRange,
     main = myTit,
     ylab = myylab,
     xlab = myxlab,
     pch=16
     )
mtext(side = 3, text = paste0(subTit))
text(x = all_ds_aucFCC,
     y = resc_auc,
     labels = names(resc_auc),
     col = curr_colors,
     cex=0.7
     )

add_curv_fit(x=all_ds_aucFCC,
             y=resc_auc, withR2 = F, lty=2)
addCorr(x=all_ds_aucFCC,
        y=resc_auc, legPos = "topright", bty="n")
legend(
  "bottomright",
  legend=names(my_colors),
  lty=-1,
  pch=16,
  col = my_colors,
  lwd = 5,
  bty="n",
  cex = 0.6)


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


############################################################################################
############################################################################################ BARPLOT
############################################################################################

plotDT <- data.frame(
  dataset = names(resc_auc_values),
  resc_AUC = unlist(resc_auc_values),
  stringsAsFactors = FALSE
)

head(plotDT)

plotDT <- plotDT[order(plotDT$resc_AUC, decreasing = TRUE),]

plotDT$dataset <- factor(plotDT$dataset, levels=as.character(plotDT$dataset))

stopifnot(plotDT$dataset %in% names(dataset_proc_colors))
curr_colors <- dataset_proc_colors[as.character(plotDT$dataset)]

myylab <- paste0(">= ", cumulVar, " thresh. resc. AUC")

mytit <- paste0(cumulVarTit, ": above thresh. resc. AUC")

mysub <- paste0("all datasets (", nrow(plotDT), ")")

p_AUC <- ggplot(plotDT, aes(x = dataset, y = resc_AUC)) +
  geom_bar(stat="identity", position="dodge", colour = barCol, fill=barCol, width=0.7) +
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(myylab),
                     breaks = scales::pretty_breaks(n = 5))+
  ggtitle(label = mytit, subtitle = mysub)+
  coord_cartesian(expand = FALSE) +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    # plot.margin = margin(t = 1, r = 1, b = 4, l = 1, unit = "pt"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(size = 6),
    axis.text.x = element_text( hjust=1,vjust = 0.5, size=10, angle = 90, color = curr_colors),
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

if(SSHFS) p_AUC

outFile <- file.path(outFold, paste0("rescAUC_", cumulVar, ".", plotType))
ggsave(p_AUC, filename = outFile, height = myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
######################################################################################
######################################################################################


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


