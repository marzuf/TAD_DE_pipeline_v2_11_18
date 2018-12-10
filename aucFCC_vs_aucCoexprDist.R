
startTime <- Sys.time()
cat(paste0("> Rscript aucFCC_vs_aucCoexprDist.R\n"))

# Rscript aucFCC_vs_aucCoexprDist.R

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("analysis_utils.R")

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


options(scipen=100)

buildTable <- T

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 5)
myWidth <- ifelse(plotType == "png", 600, 8)

myHeightScatter <- ifelse(plotType == "png", 400, 7)
myWidthScatter <- myHeightScatter

barcolors <- "gray48"

### HARD CODED
caller <- "TopDom"
script170_name <- "170_score_auc_pval_withShuffle"

famType1 <- "hgnc"
famType2 <- "hgnc_family_short"

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)


outFold <- file.path("AUCfcc_vs_AUCcoexprdist")
system(paste0("mkdir -p ", outFold))

# logFile <- file.path(outFold, paste0("auccFCC_vs_aucCoexprDist_logFile.txt"))  
# system(paste0("rm -f ", logFile))

pipOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")
all_datasets <- list.files(pipOutFold)

cat(paste0("# of datasets found: ", length(all_datasets), "\n"))

myylab <- paste0("% AUC increase")

if(buildTable) {

#cat(all_datasets[56], "\n") ; stop("--ok--\n")

all_auc <- foreach(curr_dataset = all_datasets) %dopar% {
  
  aucFCC_file <- file.path(pipOutFold, curr_dataset, script170_name, "allratio_auc_pval.Rdata")
  if(!file.exists(aucFCC_file)) cat("aucFCC_file = ", aucFCC_file, "\n")
  stopifnot(file.exists(aucFCC_file))
  
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_dataset, "auc_values.Rdata")

  if(!file.exists(aucCoexprDist_file)) cat("aucCoexprDist_file = ", aucCoexprDist_file, "\n")
  stopifnot(file.exists(aucCoexprDist_file))
  
  aucCoexprDistSameFam_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller),
                                         "AUC_COEXPRDIST_WITHFAM_SORTNODUP", paste0(curr_dataset, "_",famType1),
                                         famType2, "auc_values.Rdata")
  
  if(!file.exists(aucCoexprDistSameFam_file)) cat("aucCoexprDistSameFam_file = ", aucCoexprDistSameFam_file, "\n")
  stopifnot(file.exists(aucCoexprDistSameFam_file))
  
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  
  all_aucDistSameFam <- eval(parse(text = load(aucCoexprDistSameFam_file)))
  aucCoexprDist2 <- as.numeric(all_aucDistSameFam["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist2))
  aucCoexprDistSameFam <- as.numeric(all_aucDistSameFam["auc_ratio_sameFam_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDistSameFam))
  
  stopifnot(round(aucCoexprDist,4)==round(aucCoexprDist2,4))
  
  list(aucFCC = aucFCC, aucCoexprDist = aucCoexprDist, aucCoexprDistSameFam=aucCoexprDistSameFam)
  
}
names(all_auc) <- all_datasets

all_auc_FCC <- sapply(all_auc, function(x) x[["aucFCC"]])
all_auc_CoexprDist <- sapply(all_auc, function(x) x[["aucCoexprDist"]])
stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist))

all_auc_CoexprDistSameFam <- sapply(all_auc, function(x) x[["aucCoexprDistSameFam"]])
stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDistSameFam))


outFile <- file.path(outFold, "all_auc_FCC.Rdata")
save(all_auc_FCC, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, "all_auc_CoexprDist.Rdata")
save(all_auc_CoexprDist, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_auc_FCC.Rdata")
  all_auc_FCC <- eval(parse(text = load(outFile)))
  outFile <- file.path(outFold, "all_auc_CoexprDist.Rdata")
  all_auc_CoexprDist <- eval(parse(text = load(outFile)))
}

##############################################################################################
############################################## SCATTERPLOT ALL DATASETS - coexpr dist same Fam vs coexpr dist
##############################################################################################

stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist) )
stopifnot(names(all_auc_CoexprDistSameFam) == names(all_auc_CoexprDist) )

cat(names(all_auc_FCC)[!names(all_auc_FCC) %in% names(dataset_proc_colors) ], "\n")

stopifnot(names(all_auc_FCC) %in% names(dataset_proc_colors) )

# to plot % of AUC increase
all_auc_FCC <- (all_auc_FCC-1) * 100
all_auc_CoexprDist <- (all_auc_CoexprDist-1) * 100
all_auc_CoexprDistSameFam <- (all_auc_CoexprDistSameFam-1) * 100

# myxlab <- "AUC ratio FCC"
# myylab <- "AUC ratio coexpr."
myTit <- "AUC coexpr. same fam. vs. AUC coexpr. same fam."
myxlab <- paste0("% coexpr. AUC increase")
myylab <- paste0("% coexpr. same fam AUC increase")
myTit <- paste0("% AUC increase coexpr. same fam. vs. % AUC increase coexpr.")


curr_colors <- dataset_proc_colors[names(all_auc_CoexprDist)]

outFile <- file.path(outFold, paste0("aucCoexprDistSameFam_vs_aucCoexprDist_all_datasets.", plotType))
do.call(plotType, list(outFile, height = myHeightScatter, width =myWidthScatter))
plot(x = all_auc_CoexprDist,
     y = all_auc_CoexprDistSameFam,
#     xlim = range(all_auc_FCC),#*c(0.95, 1.15),
#     ylim = range(all_auc_CoexprDist),#*c(0.95,1.15),
     # xlim = range(all_auc_FCC)+c(-0.05,0.05),#*c(0.95, 1.15),
     # ylim = range(all_auc_CoexprDist)+c(-0.05,0.05),#*c(0.95,1.15),
    xlim = range(all_auc_CoexprDist)+c(-5,5),#*c(0.95, 1.15),
    ylim = range(all_auc_CoexprDistSameFam)+c(-5, 5),#*c(0.95,1.15),
     xlab = paste0(myxlab),
     ylab = paste0(myylab),
     pch = 16, cex=0.7,
     col = curr_colors,
     main = myTit)
mtext(text = paste0(caller, " - # of datasets = ", length(all_auc_FCC)), side = 3)
text(x = all_auc_CoexprDist,
     y = all_auc_CoexprDistSameFam,
     labels = names(all_auc_CoexprDist),
     col = curr_colors,
     pos=3, cex = 0.7)

addCorr(x=all_auc_CoexprDist, 
       y=all_auc_CoexprDistSameFam, 
       legPos="topright", 
       corMet="spearman",
       bty="n") 

# corTest <- cor.test(all_auc_FCC, all_auc_CoexprDist)
# legTxt <- paste0("PCC = ", round(corTest$estimate,2), "\n(p-val =  ", sprintf("%1.2e", corTest$p.value), ")")
#legend("topleft", legend = legTxt, bty="n")
# legend("topright", legend = legTxt, bty="n")

my_colors_leg <- my_colors
#names(my_colors_leg)[names(my_colors_leg) == "psychiatric disorder"] <- "psychiatric\ndisorder"
#names(my_colors_leg)[names(my_colors_leg) == "embryonic development"] <- "embryonic\ndevelopment"

legend("bottomright",
       legend=names(my_colors_leg),
       lty=1,
       col = my_colors_leg,
       lwd = 5,
       bty="n",
       cex = 0.7)

abline(h = 0, lty=2, col="darksalmon")
abline(v = 0, lty=2, col="darksalmon")


# abline(lm(all_auc_CoexprDist ~ all_auc_FCC), col="grey", lty=2)
#add_curv_fit(x = all_auc_CoexprDist, y=all_auc_CoexprDistSameFam, withR2 = TRUE, R2shiftX = -0.03, R2shiftY = 0, col="grey", lty=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




#stop("--ok\n")

##############################################################################################
############################################## SCATTERPLOT ALL DATASETS
##############################################################################################

stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist) )
stopifnot(names(all_auc_FCC) %in% names(dataset_proc_colors) )

# to plot % of AUC increase
all_auc_FCC <- (all_auc_FCC-1) * 100
all_auc_CoexprDist <- (all_auc_CoexprDist-1) * 100

# myxlab <- "AUC ratio FCC"
# myylab <- "AUC ratio coexpr."
myTit <- "AUC coexpr. vs. AUC FCC"
myxlab <- paste0("% FCC AUC increase")
myylab <- paste0("% coexpr. AUC increase")
myTit <- paste0("% AUC increase coexpr. vs. FCC")


curr_colors <- dataset_proc_colors[names(all_auc_FCC)]

outFile <- file.path(outFold, paste0("aucFCC_vs_aucCoexprDist_all_datasets.", plotType))
do.call(plotType, list(outFile, height = myHeightScatter, width =myWidthScatter))
plot(x = all_auc_FCC,
     y = all_auc_CoexprDist,
     #     xlim = range(all_auc_FCC),#*c(0.95, 1.15),
     #     ylim = range(all_auc_CoexprDist),#*c(0.95,1.15),
     # xlim = range(all_auc_FCC)+c(-0.05,0.05),#*c(0.95, 1.15),
     # ylim = range(all_auc_CoexprDist)+c(-0.05,0.05),#*c(0.95,1.15),
     xlim = range(all_auc_FCC)+c(-5,5),#*c(0.95, 1.15),
     ylim = range(all_auc_CoexprDist)+c(-5, 5),#*c(0.95,1.15),
     xlab = paste0(myxlab),
     ylab = paste0(myylab),
     pch = 16, cex=0.7,
     col = curr_colors,
     main = myTit)
mtext(text = paste0(caller, " - # of datasets = ", length(all_auc_FCC)), side = 3)
text(x = all_auc_FCC,
     y = all_auc_CoexprDist,
     labels = names(all_auc_FCC),
     col = curr_colors,
     pos=3, cex = 0.7)

addCorr(x=all_auc_FCC, 
        y=all_auc_CoexprDist, 
        legPos="topright", 
        corMet="spearman",
        bty="n") 

# corTest <- cor.test(all_auc_FCC, all_auc_CoexprDist)
# legTxt <- paste0("PCC = ", round(corTest$estimate,2), "\n(p-val =  ", sprintf("%1.2e", corTest$p.value), ")")
#legend("topleft", legend = legTxt, bty="n")
# legend("topright", legend = legTxt, bty="n")

my_colors_leg <- my_colors
#names(my_colors_leg)[names(my_colors_leg) == "psychiatric disorder"] <- "psychiatric\ndisorder"
#names(my_colors_leg)[names(my_colors_leg) == "embryonic development"] <- "embryonic\ndevelopment"

legend("topleft",
       legend=names(my_colors_leg),
       lty=1,
       col = my_colors_leg,
       lwd = 5,
       bty="n",
       cex = 0.7)

# abline(lm(all_auc_CoexprDist ~ all_auc_FCC), col="grey", lty=2)
add_curv_fit(x = all_auc_FCC, y=all_auc_CoexprDist, withR2 = TRUE, R2shiftX = -0.03, R2shiftY = 0, col="grey", lty=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


##############################################################################################
############################################################################################## BARPLOT WITH BOTH AUC RATIOS
##############################################################################################

stopifnot(names(all_auc_FCC) == names(all_auc_CoexprDist) )

auc_DT <- data.frame(dataset=names(all_auc_FCC),
                     auc_fcc = all_auc_FCC,
                     auc_coexpr = all_auc_CoexprDist,
                     stringsAsFactors = FALSE)
rownames(auc_DT) <- NULL

auc_DT <- auc_DT[order(auc_DT$auc_fcc, decreasing = TRUE),]

auc_DT_m <- melt(auc_DT, id=c("dataset"))
auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))

stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]

plotDT <- auc_DT_m
# plotDT$value <- plotDT$value - 1

my_breaks <- scales::pretty_breaks(n = 5)(plotDT$value)
my_labels <- my_breaks
# my_labels <- my_breaks + 1

# myylab <- paste0("AUC ratio")
myylab <- "% AUC increase"
mytit <- paste0(myylab, " - FCC and coexpr")

p_AUC <- ggplot(plotDT, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(myylab),
                     breaks = my_breaks,
                     labels = my_labels
                     # breaks = scales::pretty_breaks(n = 5)
                     ) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  ggtitle(label = mytit) +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."))+
  labs(fill  = "") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(size = 6),
    axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
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

outFile <- file.path(outFold, paste0("aucFCC_aucCoexprDist_all_datasets_barplot_orderFCC.", plotType))
ggsave(p_AUC, filename = outFile, height = 6, width=12)
cat(paste0("... written: ", outFile, "\n"))



# auc_DT <- auc_DT[order(auc_DT$auc_coexpr, decreasing = TRUE),]
# auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))

auc_DT <- auc_DT[order(auc_DT$auc_coexpr, decreasing = TRUE),]
plotDT$dataset <- factor(as.character(plotDT$dataset), levels = as.character(auc_DT$dataset))

my_breaks <- scales::pretty_breaks(n = 5)(plotDT$value)
# my_labels <- my_breaks + 1
my_labels <- my_breaks

p_AUC <- ggplot(plotDT, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle(label = mytit) +
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(myylab),
                     breaks = my_breaks,
                     labels = my_labels
                     # breaks = scales::pretty_breaks(n = 5)
                     ) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(auc_fcc = "dodgerblue4", auc_coexpr = "darkorange2"),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."))+
  labs(fill  = "") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(size = 6),
    axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
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
  )#+
  # geom_hline(yintercept = 1, linetype = 2)

if(SSHFS) p_AUC

outFile <- file.path(outFold, paste0("aucFCC_aucCoexprDist_all_datasets_barplot_orderCoexpr.", plotType))
ggsave(p_AUC, filename = outFile, height = 6, width=12)
cat(paste0("... written: ", outFile, "\n"))





##############################################################################################
############################################################################################## BARPLOT WITH ONE AUC RATIOS
##############################################################################################

auc_DT <- auc_DT[order(auc_DT$auc_fcc, decreasing = TRUE),]

auc_DT_m <- melt(auc_DT, id=c("dataset"))
auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))

stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]

plotDT <- auc_DT_m[auc_DT_m$variable == "auc_fcc",]
stopifnot(nrow(plotDT) > 0)
# plotDT$value <- plotDT$value - 1


myylab <- paste0("% AUC increase - FCC")
# myTit <- paste0("AUC ratio - FCC")
myTit <- paste0("% AUC increase - FCC")
my_breaks <- scales::pretty_breaks(n = 5)(plotDT$value)
# my_labels <- my_breaks + 1
my_labels <- my_breaks

p_AUC <- ggplot(plotDT, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat="identity", position="dodge", width = 0.7) +
  scale_x_discrete(name="")+
  ggtitle(label=myTit) +
  scale_y_continuous(name=paste0(myylab),
                     breaks = my_breaks,
                     labels = my_labels
                     # breaks = scales::pretty_breaks(n = 5)
                     ) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(auc_fcc = barcolors, auc_coexpr = barcolors),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."),
                    guide = FALSE )+
  labs(fill  = "") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
#    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(size = 6),
    axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
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

outFile <- file.path(outFold, paste0("aucFCC_all_datasets_barplot.", plotType))
ggsave(p_AUC, filename = outFile, height = myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))



# auc_DT <- auc_DT[order(auc_DT$auc_coexpr, decreasing = TRUE),]
# auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))

auc_DT <- auc_DT[order(auc_DT$auc_coexpr, decreasing = TRUE),]
auc_DT$dataset <- factor(as.character(auc_DT$dataset), levels = as.character(auc_DT$dataset))
auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT$dataset))

plotDT <- auc_DT_m[auc_DT_m$variable == "auc_coexpr",]
stopifnot(nrow(plotDT) > 0)
# plotDT$value <- plotDT$value - 1

curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]

myylab <- paste0("% AUC increase - coexpr.")

my_breaks <- scales::pretty_breaks(n = 5)(plotDT$value)
# my_labels <- my_breaks + 1
my_labels <- my_breaks

p_AUC <- ggplot(plotDT, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat="identity", position="dodge", width = 0.7) +
  ggtitle(label=paste0(myylab)) + 
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(myylab),
                     breaks = my_breaks,
                     labels = my_labels
                     # breaks = scales::pretty_breaks(n = 5)
                     ) +#, limits = c(0, max(auc_DT_m$value)+0.05))+
  # facet_grid(~dataset, switch="x") +
  coord_cartesian(expand = FALSE) +
  scale_fill_manual(values = c(auc_fcc = barcolors, auc_coexpr = barcolors),
                    labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."),
                    guide = FALSE)+
  labs(fill  = "") +
  theme( # Increase size of axis lines
    # top, right, bottom and left
#    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(size = 6),
    axis.text.x = element_text( hjust=1,vjust = 0.5, size=6, angle = 90, color = curr_colors),
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
  )#+
  # geom_hline(yintercept = 1, linetype = 2)

if(SSHFS) p_AUC

outFile <- file.path(outFold, paste0("aucCoexprDist_all_datasets_barplot.", plotType))
ggsave(p_AUC, filename = outFile, height = myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))






