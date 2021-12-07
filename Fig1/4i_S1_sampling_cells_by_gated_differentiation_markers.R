###################################################################################################
# (4i experiment) Defining single-cell differentiation states based on four differentiation markers
# Natacha Comandante-Lou
# Nov, 2021
###################################################################################################
rm(list=ls()) #clear all
cat("\014")  #clc
setwd("~/Code/Fig1")
## Libraries Required-------------------------------------------------
library(tidyverse)
library(reshape2)
library(ggpubr)
library(purrr)
library(jcolors)
library(ggsci)
library(pals)
library(ggbeeswarm)


# Read Data (Supplemental_Table_S4_4i_Baseline_Single-cell_Protein_Log_Measurements.csv) ------------
filename = '/Volumes/GoogleDrive/My\ Drive/Fallahi\ Lab/AP-1\ Project/20201002_20_Lines_Analysis/Gating_Differentiation_States/20210729_Gating_Diff_final/20201002_20Line_24h-120h_log_selected_columns_in_baseline_analysis.csv'
XXX = read.csv(filename) #data in natural log scale
XXX$condition_id = as.factor(XXX$condition_id)
XXX$timepoint_id = as.factor(XXX$timepoint_id)


## remove inf row
XXX = do.call(data.frame,lapply(XXX, function(x) replace(x, is.infinite(x),NA)))
ir = rowSums(is.na(as.matrix(XXX)))>0
X.dmso = XXX[!ir,]

# Gating differentiation markers --------
cutoffs = data.frame(AXL = log(10^2.43), NGFR =log(10^2.0), MiTF = log(10^3.2), Sox10 = log(10^2.96))

Diff = c('AXL','NGFR','MiTF','Sox10')
for (m in Diff){
  X.dmso[[sprintf('%s.class',m)]] = factor(X.dmso[[m]]>cutoffs[[m]],levels=c("TRUE","FALSE"),labels=c("H","L"))
}

## count cells from each group
df = select(X.dmso,c('MiTF.class','Sox10.class','NGFR.class','AXL.class'))
df %>% group_by(AXL.class,NGFR.class,MiTF.class,Sox10.class) %>% summarise(count=n())


# Assign differentiation class based on Tsoi et al. 2018 definition-----------
set.seed(1)
df.M = filter(X.dmso,MiTF.class%in%'H' & Sox10.class%in%'H' & NGFR.class%in%'L' & AXL.class%in%'L')%>% mutate(Diff.class = 'M') #Melanocytic
df.T = filter(X.dmso,MiTF.class%in%'H' & Sox10.class%in%'H' & NGFR.class%in%'H' & AXL.class%in%'L')%>% mutate(Diff.class = 'T') #Transitory
df.N = filter(X.dmso,MiTF.class%in%'L' & Sox10.class%in%'H' & NGFR.class%in%'H' & AXL.class%in%'H')%>% mutate(Diff.class = 'N') #Neural crest-like
df.NU = filter(X.dmso,MiTF.class%in%'L' & Sox10.class%in%'L' & NGFR.class%in%'H' & AXL.class%in%'H')%>% mutate(Diff.class = 'U.NGFR_High') #(removed from analysis since it's dominated by one cell line)
df.U = filter(X.dmso,MiTF.class%in%'L' & Sox10.class%in%'L' & NGFR.class%in%'L' & AXL.class%in%'H')%>% mutate(Diff.class = 'U.NGFR_Low') #Undifferentiated 

df = rbind(df.M,df.T,df.N,df.NU,df.U)
df$Diff.class = factor(df$Diff.class,levels = c("U.NGFR_Low", "U.NGFR_High", "N", "T", "M") )


## count cell-line by differentiation class
cellline.class.count = group_by(df,Diff.class,cellline)%>%tally()
class.count =group_by(df,Diff.class)%>%tally()

# Recursive selection of cells from different cell-lines---------------------

df.sample = filter(df,Diff.class%in%"U.NGFR_High")%>%ungroup()
for (class in c("U.NGFR_Low","N","T","M")){
  
  set.seed(1000)
  max_n= 2500 #maximum number of cells per class
  n = 0
  
  count =  filter(cellline.class.count,Diff.class%in%class)
  samples_per_group = count
  samples_per_group$n = 0
  avg_count = round((max_n-n)/ sum(count$n>0))
  cellline_sampling_func = function (){
    avg_count = round((max_n-n)/ sum(count$n>0))
    #total number of cells sampled
    samples_per_group$n[count$n>avg_count] = samples_per_group$n[count$n>avg_count] + avg_count
    samples_per_group$n[(count$n<avg_count)&(count$n>0)] = samples_per_group$n[(count$n<avg_count)&(count$n>0)] + count$n[(count$n<avg_count)&(count$n>0)]
    #number of cells left
    count$n[(count$n<avg_count)&(count$n>0)] = 0
    count$n[count$n>avg_count] = count$n[count$n>avg_count] -avg_count
    
    n = sum(samples_per_group$n)
    assign("samples_per_group", samples_per_group, envir = .GlobalEnv)
    assign("count", count, envir = .GlobalEnv)
    assign("n", n, envir = .GlobalEnv)
    assign("avg_count",avg_count,envir = .GlobalEnv)

  }
  
  
  while(sum(samples_per_group$n)<max_n&avg_count>0){
    cellline_sampling_func()
  }
  
  
  dff = filter(df,Diff.class%in%class)%>%group_by(cellline)%>%slice(sample(n(), min(max(samples_per_group$n), n()))) %>% ungroup()
  df.sample = rbind(df.sample,dff)
}

# Plot sampled cells ------------------------------------
pal1 = pal_d3(alpha = 0.55)(8)
pal2 = pal_rickandmorty(palette = c("schwifty"), alpha = 0.9)(12)
mypal = c(pal1,pal2)

ggplot(df.sample)+
  geom_bar(aes(y = Diff.class,fill = cellline))+
  guides(colour = guide_legend(override.aes = list(size=5),title = "Cell Lines" ))+
  theme_classic() +
  scale_fill_manual(values = mypal)+
  ylab("Differentiaion Class") +
  ggtitle("Cell-line Distributions Across Differentiation States")

ggsave("cell-line_distributation_across_diff_states.pdf", width = 8.2, height = 5.5)


# Save sampled cells which will be used for further analysis
df.sample = filter(df.sample,!Diff.class%in%c('U.NGFR_High')) #(removed U.NGFR_high from analysis since it's dominated by one cell line)
write_csv(df.sample,'Sampled_Gated_Cells_by_Diff_States.csv') #Used for random forest model

