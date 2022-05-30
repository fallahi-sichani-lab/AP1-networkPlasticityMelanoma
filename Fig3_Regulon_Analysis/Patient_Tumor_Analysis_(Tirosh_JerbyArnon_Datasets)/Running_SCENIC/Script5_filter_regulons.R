###################################################################################################
# (SCENIC) Filter and summarize regulons from 100x SCENIC runs
# Adapted from Wouters et al 2020
# https://github.com/aertslab/singlecellRNA_melanoma_paper
# Nov, 2021
###################################################################################################
setwd("~/Code/Fig3/")
library(data.table)
library(AUCell)
library(hdf5r)

nRuns = 100
regulons <- list()


setwd(staging)

# [1] get regulons from 100x runs
###################################### {
mtfdir <- file.path(staging, sprintf("scenic-output-%dx",nRuns))

# 100x pyscenic runs
run <- "mlt"
dir <- "pos"

db <- "mtf"
for (i in 1:nRuns) {
  reg <- readRDS(file.path(mtfdir, paste0("run_", i), paste0("reg_", db, ".rds.gz")))
  regulons[[paste(run, i, db, dir, sep = "_")]] <- reg
}

#remove empty lists
regulons <- regulons[!lapply(regulons, length)==0] 


saveRDS(regulons, file = "post-analysis/regulons_extended_mlt.rds.gz", compress = "gzip")

# [2]  make one data frame with presence/absence
################################################ {
reg_long <- lapply(regulons, function(regset) {
stack(regset$pos, drop = TRUE)
})
reg_long <- rbindlist(reg_long, idcol = names(reg_long))
colnames(reg_long) <- c("psrun", "target", "tf")
reg_long[, regulon := paste0(psrun, "_", tf)]
incMat <- as.data.frame(with(reg_long, table(target, regulon)) > 0L) + 0L
dim(incMat)
#}

# [3] summarize all runs
######################## {
## 100x scenic runs

#how often is tf-target connection found in 100x runs
mlt_runs_mtf <- sapply(1:nRuns, function(i) {paste0("mlt_", i, "_mtf_pos")})
#how often is regulator found in 100x runs
tf_recurrence <- unique(reg_long[psrun %in% mlt_runs_mtf, .(psrun, tf)])[, .N, by = tf]
colnames(tf_recurrence) <- c("tf", "tf_rec_100x_mtf")

#join both
mlt_runs_mtf <- reg_long[psrun %in% mlt_runs_mtf, .N, by = .(target, tf)]
colnames(mlt_runs_mtf) <- c("target", "tf", "conn_rec_100x_mtf")

mlt_runs_mtf[, dir := "pos"]
mlt_runs_mtf <- mlt_runs_mtf[tf_recurrence, on = "tf"]


#full outer join

allruns = mlt_runs_mtf


saveRDS(allruns, file = "post-analysis/all_runs_extended_summary_mtf.rds.gz", compress = "gzip")



# [4] filter regulons & calculate AUCell
########################################################################### {

allruns <- allruns[, .(tf, target, dir,
                       conn_rec_100x_mtf, tf_rec_100x_mtf)]
for (i in names(allruns)) allruns[is.na(get(i)), (i):=0]


#filter rule:
#keep only positive TF-target connections
#for 100x recurrent TFs, take targets that come up 80x
#for at least 80x recurrent TFs, take all targets


filt_mtf <- with(allruns,
	dir == "pos" &
	#!tf %in% nontfs &
	((tf_rec_100x_mtf == nRuns & conn_rec_100x_mtf >= 80) |
	(tf_rec_100x_mtf < nRuns & tf_rec_100x_mtf >= 80 & conn_rec_100x_mtf > 0)))


regs_mtf <- sapply(as.vector(unique(allruns[filt_mtf, tf])), function(TF) {
    	as.vector(allruns[filt_mtf & tf == TF, target])
}, simplify = FALSE, USE.NAMES = TRUE)
names(regs_mtf) <- paste0(names(regs_mtf), "_mtf")

regs <- c(regs_mtf)

summary(unlist(lapply(regs, length)))

# rm regulons <10 targets
regs <- regs[!unlist(lapply(regs, length)) < 10]

# calculate AUCell
#get gene expression matrix from loom file 
loom <- H5File$new("pre-processing/jerby-arnon_sc_skin_filtered_for_scenic.loom", mode = "r") 
dgem <- t(loom[["matrix"]]$read()) 
genes <- loom[["row_attrs"]][["Gene"]]$read() 
rownames(dgem) <- genes 
cells <- loom[["col_attrs"]][["CellID"]]$read() 
colnames(dgem) <- cells 
loom$close_all() 

#calculate AUCell 
#use ranking from make_AUCell_ranking_db.R
aucellRankings <- readRDS("aucellRankings.rds.gz")
regulonAUC <- AUCell_calcAUC(regs, aucellRankings,  
    	aucMaxRank = aucellRankings@nGenesDetected["1%"], nCores = 1) 
auc <- getAUC(regulonAUC) 

saveRDS(regulonAUC, file = "post-analysis/all_pos_regulons_recurrent_100x_mtf_regulonAUC01.rds.gz", compress = "gzip")
saveRDS(auc, file = "post-analysis/all_pos_regulons_recurrent_100x_mtf_aucell01.rds.gz", compress = "gzip") 
saveRDS(regs, file = "post-analysis/all_pos_regulons_recurrent_100x_mtf_regulons.rds.gz", compress = "gzip") 

#save as text file
rownames(auc) <- gsub("_mtf", "", rownames(auc))

auc <- as.data.table(auc, keep.rownames = "rn")
fwrite(auc, file = "post-analysis/all_pos_regulons_recurrent_100x_mtf_aucell01.tsv", quote = F, sep = "\t")
#}


sessionInfo() 
