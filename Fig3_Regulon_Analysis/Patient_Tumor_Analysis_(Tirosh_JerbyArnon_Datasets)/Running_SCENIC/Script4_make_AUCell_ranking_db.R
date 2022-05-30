###################################################################################################
# (SCENIC) Make data base with AUCell rankings for Tirosh et al Jerby-Arnon data
# Adapted from Wouters et al 2020
# https://github.com/aertslab/singlecellRNA_melanoma_paper
# Nov, 2021
###################################################################################################
setwd("~/Code/Fig3/")
library(hdf5r)
library(AUCell)
seed <- 123


loomfile <- "pre-processing/jerby-arnon_sc_skin_filtered_for_scenic.loom"
loom <- H5File$new(loomfile, mode = "r")

### retrieve gene expression matrix (genes x cells)
dgem <- t(loom[["matrix"]]$read())

### get genes and cells
genes <- loom[["row_attrs"]][["Gene"]]$read()
rownames(dgem) <- genes
cells <- loom[["col_attrs"]][["CellID"]]$read()
colnames(dgem) <- cells

### close loom file
loom$close_all()

# create rankings 
set.seed(seed)
aucellRankings <- AUCell_buildRankings(dgem, nCores = 20, plotStats = F)
saveRDS(aucellRankings, file = "aucellRankings.rds.gz", compress = "gzip")


sessionInfo() 