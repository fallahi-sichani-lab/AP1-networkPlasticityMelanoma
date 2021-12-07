###################################################################################################
# (SCENIC) Calculate enrichment of Tsoi differentiation signatures (AUCell) for 
# individual cells from Tirosh et al and Jerby-Arnon et al scRNA-Seq data
# Natacha Comandante-Lou
# Nov, 2021
###################################################################################################
rm(list=ls()) #clear all
cat("\014")  #clc
setwd("~/Code/Fig3")

#Vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html#some-tips
library(AUCell)
library(AUCell)
library(SCopeLoomR)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(colorspace)
library(jcolors)
library(ggsci)
library(pals)
library(RColorBrewer)
library(gplots)
###################################################
# Get exprMat from loom file --------------------------
# jerby_tirosh et al

jerby_tirosh_pyScenicLoomFile <- "jerby-arnon_sc_skin_filtered_for_scenic.loom"
jerby_tirosh_loom <- open_loom(jerby_tirosh_pyScenicLoomFile, mode="r")

jerby_tirosh_exprMat <- get_dgem(jerby_tirosh_loom)
jerby_tirosh_cellInfo <- get_cell_annotation(jerby_tirosh_loom)
jerby_tirosh_df.metadata = data.frame(cell_id = rownames(jerby_tirosh_cellInfo),tumor = jerby_tirosh_cellInfo$Tumor)


# 1. Define Gene Sets ----------------------------------------------------------------------------------------

# Overlapping gene_set-------------
grouping = "M-MT_NU-U"
M_MT.genes <- scan("Tsoi_Gene_Set/M_MT_Gene_list.txt", what="", sep="\n")

NU_U.genes <- scan("Tsoi_Gene_Set/NU_U_Gene_list.txt", what="", sep="\n")

DiffGeneSets <- c(
  GeneSet(M_MT.genes, setName="M_MT"),
  GeneSet(NU_U.genes, setName="NU_U")
  
)


geneSets <- GeneSetCollection(DiffGeneSets)
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep="")) #add the gene-set size into its name

# 2. Build gene-expression rankings for each cell ---------------------------------------------------------
#For each cell, the genes are ranked from highest to lowest value. The genes with same expression value are shuffled.

#check that most cells have at least the number of 
#expressed/detected genes that are going to be used to calculate the AUC
par(mfrow=c(3,1)) 

cells_rankings <- AUCell_buildRankings(jerby_tirosh_exprMat, nCores=1, plotStats=TRUE)
#Quantiles for the number of genes detected by cell: 
# (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min      1%      5%     10%     50%    100% 
# 1933.0  2237.4  2938.0  3511.6  5444.0 11106.0 


# 3. Calculate enrichment for the gene signatures (AUC) -----------------------------------------------------
#To determine whether the gene set is enriched at the top of the gene-ranking for each cell, 
#AUCell uses the “Area Under the Curve” (AUC) of the recovery curve.

#Tsoi_diff_cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
Tsoi_diff_cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = cells_rankings@nGenesDetected[["1%"]])
# Genes in the gene sets NOT available in the dataset: 
# M_MT (187g): 	8 (4% of 187)
# NU_U (224g): 	103 (46% of 224)


# create dataframe for results
dff =  data.frame(t(Tsoi_diff_cells_AUC@assays@data@listData[["AUC"]]))
colnames(dff) = Tsoi_diff_cells_AUC@NAMES
dff$cell_id = rownames(dff)

df.tsoi_diff_auc = plyr::join(jerby_tirosh_df.metadata,dff,by = "cell_id")


# 4. Determine the cells with the given gene signatures or active gene sets---------------------------------
#The thicker vertical line indicates the threshold selected by default ($aucThr$selected): 
#the highest value to reduce the false positives.
set.seed(123)

pdf(file=sprintf("jerby_tirosh_Diff_AUC_Automatic_Thresholds_%s.pdf",grouping),width = 10, height = 5)
par(mfrow=c(2,4)) 
Tsoi_diff_cells_assignment <- AUCell_exploreThresholds(Tsoi_diff_cells_AUC, plotHist=TRUE, assign=TRUE) 

dev.off()

cellsAssigned <- lapply(Tsoi_diff_cells_assignment, function(x) x$assignment)


write.csv(df.tsoi_diff_auc,sprintf("jerby_tirosh_Tsoi_Diff_AUCell_DF_%s.csv",grouping), row.names = FALSE)



