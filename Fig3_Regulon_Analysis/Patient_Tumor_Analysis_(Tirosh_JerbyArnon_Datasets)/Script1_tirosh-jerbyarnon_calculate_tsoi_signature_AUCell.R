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
grouping = "M-MT_MT-T-TN_TN-N-NU_NU-U"

M_MT.genes <- scan("Tsoi_Gene_Set/M_MT_Gene_list.txt", what="", sep="\n")
MT_T_TN.genes <-scan("Tsoi_Gene_Set/MT_T_TN_Gene_list.txt", what="", sep="\n")
TN_N_NU.genes <-scan("Tsoi_Gene_Set/TN_N_NU_Gene_list.txt", what="", sep="\n")
NU_U.genes <- scan("Tsoi_Gene_Set/NU_U_Gene_list.txt", what="", sep="\n")

DiffGeneSets <- c(
  GeneSet(M_MT.genes, setName="M_MT"),
  GeneSet(MT_T_TN.genes, setName="MT_T_TN"),
  GeneSet(TN_N_NU.genes, setName="TN_N_NU"),
  GeneSet(NU_U.genes, setName="NU_U")
  
)



#--------------------------------------
geneSets <- GeneSetCollection(DiffGeneSets)
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep="")) #add the gene-set size into its name

# 2. Build gene-expression rankings for each cell ---------------------------------------------------------
#For each cell, the genes are ranked from highest to lowest value. The genes with same expression value are shuffled.

#check that most cells have at least the number of 
#expressed/detected genes that are going to be used to calculate the AUC
par(mfrow=c(3,1)) 

cells_rankings <- AUCell_buildRankings(jerby_tirosh_exprMat, nCores=1, plotStats=TRUE)
# Quantiles for the number of genes detected by cell: 
#   (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min       1%       5%      10%      50%     100% 
# 1703.00  1883.62  2274.55  2747.10  5315.00 12471.00 


# 3. Calculate enrichment for the gene signatures (AUC) -----------------------------------------------------
#To determine whether the gene set is enriched at the top of the gene-ranking for each cell, 
#AUCell uses the “Area Under the Curve” (AUC) of the recovery curve.

#Tsoi_diff_cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
Tsoi_diff_cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = cells_rankings@nGenesDetected[["1%"]])
# Genes in the gene sets NOT available in the dataset: 
#   M_MT (187g): 	9 (5% of 187)
# MT_T_TN (179g): 	14 (8% of 179)
# TN_N_NU (197g): 	58 (29% of 197)
# NU_U (224g): 	96 (43% of 224)


# create dataframe for results
dff =  data.frame(t(Tsoi_diff_cells_AUC@assays@data@listData[["AUC"]]))
colnames(dff) = Tsoi_diff_cells_AUC@NAMES
dff$cell_id = rownames(dff)

df.tsoi_diff_auc = plyr::join(jerby_tirosh_df.metadata,dff,by = "cell_id")


# 4. Determine the cells with the given gene signatures or active gene sets---------------------------------
#The thicker vertical line indicates the threshold selected by default ($aucThr$selected): 
#the highest value to reduce the false positives.


pdf(file=sprintf("jerby_tirosh_Diff_AUC_Automatic_Thresholds_%s.pdf",grouping),width = 10, height = 5)
par(mfrow=c(2,4)) 
Tsoi_diff_cells_assignment <- AUCell_exploreThresholds(Tsoi_diff_cells_AUC, plotHist=TRUE, assign=TRUE) 

dev.off()

cellsAssigned <- lapply(Tsoi_diff_cells_assignment, function(x) x$assignment)


write.csv(df.tsoi_diff_auc,sprintf("jerby_tirosh_Tsoi_Diff_AUCell_DF_%s.csv",grouping), row.names = FALSE)



