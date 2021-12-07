###################################################################################################
# (SCENIC) Get regulons from pyscenic output table and top motifs for each regulon (max NES)
# Adapted from Wouters et al 2020
# 
# Nov, 2021
###################################################################################################
setwd("~/Code/Fig3/")

library(data.table) 
library(parallel)


scenicdir <- file.path('../scenic-output-100x')
setwd(scenicdir)
nRuns = 100
# function to get all target genes for a specific TF from reg table  
getRegulon <- function(reg, tf, d) {  
   #d: "neg"/"pos" 
   digitpattern <- "[[:digit:]]+\\.[[:digit:]e-]+"
   alltargets <- lapply(1:reg[TF == tf & dir == d, .N], function(i) {  
        targets <- strsplit(  
         gsub(digitpattern, "",  
          gsub(" ", "",  
           gsub("'", "",   
            gsub("]", "",  
             gsub("[", "",   
              gsub("[()]", "", reg[TF == tf & dir == d, TargetGenes][i]),   
             fixed = TRUE),   
            fixed = TRUE),  
           )  
          )  
         )  
        , ",")[[1]]  
   targets[sapply(targets, function(p) {p != ""})]  
   })  
   Reduce(union, alltargets)  
}

mclapply(1:nRuns, function(runno) {
	motiffile <- paste0("run_", runno, "/","motif_reg_nodom.csv")
        tab <- fread(motiffile, skip = 3, sep = ",")
        colnames(tab) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity",
			                  "Annotation", "Context", "TargetGenes", "RankAtMax")

	tab[, dir := ifelse(grepl("activating", Context), "pos", "neg")]

	regulons <- list()
	d <- "neg"
	regulons[[d]] <- sapply(unique(tab[dir == d, TF]), function(tf) {  
		getRegulon(tab, tf, d)})  
	d <- "pos"
	regulons[[d]] <- sapply(unique(tab[dir == d, TF]), function(tf) {  
		getRegulon(tab, tf, d)})
	
	print("final regulons:")
	print(sapply(regulons, length))
	
	saveRDS(regulons, file = paste0("run_", runno, "/reg_mtf.rds.gz"), compress = "gzip")

        tab <- tab[, .(TF, MotifID, NES)]
        tab <- tab[tab[, .I[NES == max(NES)], by = TF]$V1]

        saveRDS(tab, file = paste0("run_", runno, "/top_mtf_per_TF.rds.gz"), compress = "gzip")
}, mc.cores = 30)

#combine top motifs into one table
motifs <- rbindlist(mclapply(1:nRuns, function(runno) {
	tab <- readRDS(file.path(paste0("run_", runno, "/top_mtf_per_TF.rds.gz")))
	tab[, run := runno]
	tab
}, mc.cores = 30))
saveRDS(motifs, file = sprintf("top_mtf_per_TF_%druns.rds.gz",nRuns), compress = "gzip")

sessionInfo() 

