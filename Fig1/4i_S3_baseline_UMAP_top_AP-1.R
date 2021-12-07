###################################################################################################
# (4i experiment) UMAP based on top AP-1 factors
# Natacha Comandante-Lou
# Nov, 2021
###################################################################################################
rm(list=ls()) #clear all
cat("\014") #clc
setwd("~/Code/Fig1/")

## Libraries Required-------------------------------------------------

library(ggplot2)
library(ggpubr)
library(purrr)
library(factoextra)
library(tidyverse)
library(jcolors)
library(umap)


# List of Markers-----------------------------------------------------
AP1s = c("cFos","Phospho_cFos","Fra1","Phospho_Fra1","Fra2",
         "cJun","Phospho_cJun","JunB", "JunD",
         "Phospho_ATF1","ATF2","Phospho_ATF2","ATF3","ATF4","Phospho_ATF4","ATF5","ATF6")

selected_AP1s = c('Phospho_Fra1', 'ATF4', 'Fra2', 'Phospho_cFos', 'cFos', 'cJun') #Top6

selected_markers = c(selected_AP1s)

cc = 0

# use for loops to explore different UMAP parameters
for (nn in c(90)){ #nearest neighbors
  for (md in c(0.7)){ #minimum distance
    for (npc in c(4)){ #number of principal components
      for (metric in c("euclidean")){ #distance metrics
        # Parameters --------------------------------------------------------------
        
        set.seed(100)
        exp_name = '20Line_Baseline_Top6_AP1'
        
        ncomp = 2; #number of UMAP components
        scale_or_not = TRUE; #z-score
        seed = 42
        
        
        # Output paths ------------------------------------------------------------
        parameter = sprintf("scale=%s-%s-npc=%0.0f-nn=%0.0f-md=%0.0e-ncomp=%0.0f",
                            scale_or_not,metric,npc,nn,md,ncomp)
        print(parameter)
        output_data_path = sprintf("%s/umap_output/%s/%s",getwd(),exp_name,parameter)
        data_filename = sprintf("%s/%s.RData",output_data_path,parameter)
        
        dir.create(sprintf("%s/umap_output",getwd()))
        dir.create(sprintf("%s/umap_output/%s",getwd(),exp_name))
        dir.create(output_data_path)
        
        figure_path = sprintf("%s/umap_figure/%s/",getwd(),exp_name)
        dir.create(sprintf("%s/umap_figure",getwd()))
        dir.create(sprintf("%s/umap_figure/%s",getwd(),exp_name))
        dir.create(figure_path)
        
        umap_params = sprintf("exp_name=%s, scale=%s, metric=%s, npc=%0.0f, nn=%0.0f, md=%0.0e, ncomp=%0.0f", exp_name,scale_or_not,metric,npc,nn,md,ncomp)
        
        # Read Data --------------------------------------------------------------
       
        df.sample = read.csv('Sampled_Gated_Cells_by_Diff_States.csv')
        df.sample$Diff.class = factor(df.sample$Diff.class, levels = c("M","T","N","U.NGFR_Low"), labels = c("Melanocytic","Transitory","Neural crest-like","Undifferentiated"))
        
        #z-score
        num_col = colnames(df.sample)[sapply(df.sample, class)=="numeric"]
        df.sample.z = as.data.frame(df.sample)
        df.sample.z[,num_col] <- lapply(df.sample.z[,num_col], function(x) c(scale(x,center=TRUE,scale=scale_or_not))) 
        df.sample.z$Diff.class = droplevels(df.sample.z$Diff.class)
        # nvar = dim(df.sample)[2]
        # df.sample.z = as.data.frame(df.sample)
        # df.sample.z[,c(8:nvar-5)] <- lapply(df.sample.z[,c(8:nvar-5)], function(x) c(scale(x,center=TRUE,scale=scale_or_not))) 
        # df.sample.z$Diff.class = droplevels(df.sample.z$Diff.class)
        
        X_selected_markers <-dplyr::select(df.sample.z,all_of(selected_markers))

        # Set up R UMAP-----------------------------------------------------------
        custom.config = umap.defaults
        custom.config$random_state = seed
        custom.config$n_neighbors = as.integer(nn)
        custom.config$n_components = as.integer(ncomp)
        custom.config$metric = metric
        custom.config$min_dist = md
        print(custom.config)
        # pca
        ap1.pca = prcomp(X_selected_markers)
        ce = ap1.pca$x
        # run umap
        set.seed(seed)
        
        umap_output <- umap(as.matrix(x = ce[,1:npc]),config = custom.config)
        
        #Collect Output-------------------------------------------------------
        
        output = list(umap = umap_output$layout,
                      data = df.sample,
                      data.zscore = df.sample.z,
                      selected_marker = selected_markers,
                      parameters = data.frame(npc=npc,nn=nn,md=md,ncomp=ncomp,metric=metric,scale_or_not=scale_or_not),
                      parameter_label = umap_params)
        savefile = sprintf("%s/scale=%s-%s-npc=%0.0f-nn=%0.0f-md=%0.0e-ncomp=%0.0f.RData",output_data_path,scale_or_not,metric,npc,nn,md,ncomp)
        save(output,file = savefile)
        
        load(savefile)
        df = output[["data.zscore"]]
        df$X1=output[["umap"]][,1]
        df$X2=output[["umap"]][,2]
        Diff = c("MiTF","Sox10","NGFR","AXL")
        
        #Plot Func
        plot_embedding_categorical <- function(category,.df){
          dpal_func = jcolors_contin("pal4",reverse = TRUE,bias = 0.5)
          dpal = rev(dpal_func(length(unique(.df[[category]]))+1))
          ggplot(data = .df,aes(x = `X1`,y = `X2`,fill = .data[[category]]))+
            geom_point(color ="black",shape = 21,stroke = 0.2)+
            guides(fill= guide_legend(override.aes = list(size=3),title = "Differentiation State"))+
            scale_fill_manual(values=dpal[sort(as.numeric(unique(.df[[category]])))] )+
            xlab("UMAP 1") +
            ylab("UMAP 2") +
            ggtitle(sprintf("UMAP based on Top %d AP-1",length(selected_markers)))+
            labs(caption =parameter)+
            theme_classic()
        }
        
        
        if (cc == 0) {
          f = lapply('Diff.class',plot_embedding_categorical,df)
        } else {
          f = c(f,lapply('Diff.class',plot_embedding_categorical,df))
          
        }
        ggarrange(plotlist = f, nrow= ceiling(sqrt(length(f))), ncol = ceiling(sqrt(length(f))))
        cc = 1
        
      }
    }
  }
}

ggsave(sprintf('%s/4i_UMAP_Top_AP-1_By_Diff_Class.eps',figure_path),width = 5, height =9,limitsize = FALSE,dpi=300,device = 'eps')

