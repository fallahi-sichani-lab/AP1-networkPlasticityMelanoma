###################################################################################################
# (SCENIC) Calculate enrichment of Tsoi differentiation signatures (AUCell) for 
# individual cells from Wouters et al scRNA-Seq data
# Natacha Comandante-Lou
# Nov, 2021
###################################################################################################
rm(list=ls()) #clear all
cat("\014")  #clc
setwd("~/Code/Fig3")

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

#Vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html#some-tips

# Get exprMat from loom file ---------------------------------------------------
# Wouters et al data
# wouters_pyScenicLoomFile <- "10_Baselines_filteredRegulons.loom"
# wouters_loom <- open_loom(wouters_pyScenicLoomFile, mode="r")


wouters_loom_folder ="/Volumes/GoogleDrive/My Drive/Fallahi Lab/AP-1 Project/Wouter et al Data Analysis/Wouters_SCENIC_Output/Regulon_Regression_Model"
wouters_loom_filename = "10_Baselines_filteredRegulons.loom"
wouters_pyScenicLoomFile <- file.path(wouters_loom_folder, wouters_loom_filename)
wouters_loom <- open_loom(wouters_pyScenicLoomFile, mode="r")

wouters_exprMat <- get_dgem(wouters_loom)
wouters_cellInfo <- get_cell_annotation(wouters_loom)
wouters_df.metadata = data.frame(cell_id = rownames(wouters_cellInfo),cell_line = wouters_cellInfo$cell_line)

# 1. Define Gene Sets ----------------------------------------------------------

# Overlapping gene_set-----------------
grouping = "M-MT_MT-T-TN_NU-U"
M_MT.genes <- scan("Tsoi_Gene_Set/M_MT_Gene_list.txt", what="", sep="\n")
MT_T_TN.genes <-scan("Tsoi_Gene_Set/MT_T_TN_Gene_list.txt", what="", sep="\n")
NU_U.genes <- scan("Tsoi_Gene_Set/NU_U_Gene_list.txt", what="", sep="\n")

DiffGeneSets <- c(
  GeneSet(M_MT.genes, setName="M_MT"), #Melanocitic gene set + Melanocytic-Transitory gene set
  GeneSet(MT_T_TN.genes, setName="MT_T_TN"), #Melanocytic-Transitory gene set + Transitory gene set + Transitory-Neural crest-like gene set
  GeneSet(NU_U.genes, setName="NU_U") #Neural crest-like gene set + Undifferentiated gene set
  
)

geneSets <- GeneSetCollection(DiffGeneSets)
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep="")) #add the gene-set size into its name

# 2. Build gene-expression rankings for each cell ------------------------------
#For each cell, the genes are ranked from highest to lowest value. The genes with same expression value are shuffled.

#check that most cells have at least the number of 
#expressed/detected genes that are going to be used to calculate the AUC
par(mfrow=c(3,1)) 

cells_rankings <- AUCell_buildRankings(wouters_exprMat, nCores=1, plotStats=TRUE)
#Quantiles for the number of genes detected by cell: 
# (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min      1%      5%     10%     50%    100% 
# 1331.00 2037.31 2628.00 2837.20 3688.50 7933.00 


# 3. Calculate enrichment for the gene signatures (AUC) -------------------------
#To determine whether the gene set is enriched at the top of the gene-ranking for each cell, 
#AUCell uses the “Area Under the Curve” (AUC) of the recovery curve.

Tsoi_diff_cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = cells_rankings@nGenesDetected[["1%"]])
# Genes in the gene sets NOT available in the dataset: 
# M_MT (187g): 	5 (3% of 187)


# create dataframe for results
dff =  data.frame(t(Tsoi_diff_cells_AUC@assays@data@listData[["AUC"]]))
colnames(dff) = Tsoi_diff_cells_AUC@NAMES
dff$cell_id = rownames(dff)

df.tsoi_diff_auc = plyr::join(wouters_df.metadata,dff,by = "cell_id")


# 4. Determine the cells with the given gene signatures or active gene sets------

set.seed(123)

pdf(file=sprintf("Wouters_Diff_AUC_Automatic_Thresholds_%s.pdf",grouping),width = 10, height = 5)
par(mfrow=c(2,4)) 
Tsoi_diff_cells_assignment <- AUCell_exploreThresholds(Tsoi_diff_cells_AUC, plotHist=TRUE, assign=TRUE) 

dev.off()



write.csv(df.tsoi_diff_auc,sprintf("Wouters_Tsoi_Diff_AUCell_DF_%s.csv",grouping), row.names = FALSE)




















