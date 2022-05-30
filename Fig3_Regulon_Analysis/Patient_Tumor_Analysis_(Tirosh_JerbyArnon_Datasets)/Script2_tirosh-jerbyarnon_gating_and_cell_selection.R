###################################################################################################
# (SCENIC) Gating signature enrichment AUCell Scores to determine phenotypes of 
# individual cells from Tirosh & Jerby-Arnon et al scRNA-Seq data from patient tumors
# Natacha Comandante-Lou
# Nov, 2021
################################################################################################### 

rm(list=ls()) #clear all
cat("\014")  #clc
setwd("/~/Code/Fig3")
library(stringr)
library(ggpubr)
library(ggbeeswarm)
library(ggsci)
library(dplyr)
library(rstatix)
# Read Data --------------------------------------------------------------
grouping = "M-MT_MT-T-TN_TN-N-NU_NU-U"
X = read.csv(sprintf("JerbyArnon_Tirosh_SCENIC_Merged_Diff_Reg_AUCell_DF_%s.csv",grouping))
colnames(X) = str_replace(colnames(X),'.x','')%>%str_replace('.y','-1')%>%str_replace('.3','-1')
colnames(X) = str_replace(colnames(X),'\\.\\.',' (')%>%str_replace('\\.',')')
X$tumor = as.factor(X$tumor)

regulon_list = c("FOS_mtf","FOSB_mtf","FOSL1_mtf","FOSL2_mtf",
                 "ATF1_mtf","ATF2_mtf","ATF4_mtf", 
                 "ATF6_mtf","ATF6B_mtf","ATF7_mtf",
                 "JUN_mtf","JUNB_mtf","JUND_mtf")

sample_info = read.csv('/Volumes/GoogleDrive/My Drive/Fallahi Lab/AP-1 Project/Tirosh et al Data Analysis/Data/GSE115978/Jerby_Tirosh_Sample_Characteristics.csv')
no_pre_treatment.tumor_list = filter(sample_info,Pre.Treatment%in%"None")[["tumor"]]

x = select(X,c("cell_id","tumor","M_MT (187g)","MT_T_TN (179g)","TN_N_NU (197g)","NU_U (224g)",regulon_list))%>%filter(tumor%in%no_pre_treatment.tumor_list)
x = mutate(x, FOS_JUN_mtf_ratio = FOS_mtf/JUN_mtf )

#thresholds for each differentiation state
#1 - M_MT (187g)
#2 - MT_T_TN (179g)
#3 - TN_N_NU (197g)
#4 - NU_U (224g)

#melanocytic thresholds
thres.m.1 = 0.12
thres.m.3 = 0.04
thres.m.4 = 0.03


#transitory thresholds
thres.t.2 = 0.12
thres.t.4 = 0.04

#neural crest-like thresholds
thres.n.1 = 0.03
thres.n.3 = 0.066

#undifferentiated threshold
thres.u.1 = 0.04
thres.u.2 = 0.05
thres.u.4 = 0.045
# 



thres = data.frame(`M_MT (187g)` = c(thres.m.1,NA, thres.n.1, thres.u.1),
                  `MT_T_TN (179g)` = c(NA, thres.t.2, NA, thres.u.2),
                  `TN_N_NU (197g)` = c(thres.m.3, NA, thres.n.3, NA),
                  `NU_U (224g)` = c(thres.m.4, thres.t.4, NA, thres.u.4))
colnames(thres) = str_replace(colnames(thres),'\\.\\.',' (')%>%str_replace('\\.',')')
thres$state = c("Melanocytic","Transitory","Neural crest-like","Undifferentiated")




#Determining cell state based on their signature enrichment AUCell scores

x = mutate(x, Diff_class = case_when((x[["M_MT (187g)"]]>=thres.m.1) & (x[["TN_N_NU (197g)"]]<thres.m.3) & (x[["NU_U (224g)"]]<thres.m.4) ~ "Melanocytic",
                                     (x[["M_MT (187g)"]]<thres.m.1) & (x[["MT_T_TN (179g)"]]>=thres.t.2) & (x[["TN_N_NU (197g)"]]<thres.n.3) & (x[["NU_U (224g)"]]<thres.t.4) ~ "Transitory",
                                     (x[["M_MT (187g)"]]<thres.n.1) & (x[["MT_T_TN (179g)"]]<thres.t.2) & (x[["TN_N_NU (197g)"]]>=thres.n.3) & (x[["NU_U (224g)"]]<thres.u.4)~ "Neural crest-like",
                                     (x[["M_MT (187g)"]]<thres.u.1) & (x[["MT_T_TN (179g)"]]<thres.u.2) & (x[["NU_U (224g)"]]>=thres.u.4) ~ "Undifferentiated",
                                     TRUE ~ "NA"))



count.diff = group_by(x,Diff_class)%>%tally()
count.diff.cellline = group_by(x,Diff_class, tumor)%>%tally()

x$Diff_class = factor(x$Diff_class,levels = rev(c("Melanocytic","Transitory","Neural crest-like","Undifferentiated","NA")))

pal1 = pal_d3(palette= "category20")(12)
pal2 = pal_rickandmorty(palette = c("schwifty"), alpha = 0.9)(12)
mypal = c(pal1,pal2)



# Box plots (t-test) --------------------

.df = filter(x,!Diff_class%in%c("Transitory","Neural crest-like","NA"))%>%group_by(Diff_class) #remove T, N due to small number of cells
.df$Diff_class = factor(.df$Diff_class, levels = c("Melanocytic","Undifferentiated"))

#Create t-test comparison combinations
class_list = as.character(unique(.df$Diff_class))

combns = combn(class_list,2,simplify = FALSE)
my_comparisons <- combns

reg = regulon_list[1]

reg_boxplots_func = function(reg,.df){
  stat.test <- .df%>% ungroup()%>%t_test(comparisons = my_comparisons,
                                         formula = as.formula(sprintf("%s ~ Diff_class",reg)))%>%
    add_xy_position(x = "Diff_class", dodge = 0.9)
  
  
  ggplot(.df,aes(x = Diff_class, y = .data[[reg]],group = Diff_class)) +
    geom_quasirandom(groupOnX=TRUE,size = 1,aes(color = tumor))+
    geom_boxplot(alpha = 0.3,outlier.shape = NA)+
    scale_color_manual(values = mypal[sort(as.numeric(unique(.df[["tumor"]])),decreasing = TRUE)])+
    guides(color = guide_legend(override.aes = list(size=10),title = "tumor"))+
    #stat_compare_means(comparisons = my_comparisons,method = "t.test")+ 
    stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01)+
    ylab('AUCell')+
    
    ggtitle(reg)+
    theme(panel.background = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
}

figs = lapply(c("FOS_mtf","FOSL1_mtf","FOSL2_mtf","JUN_mtf","FOS_JUN_mtf_ratio"),reg_boxplots_func,.df)

ggarrange(plotlist = figs, nrow = 1, ncol = 5, common.legend = TRUE)
ggsave(sprintf("JerbyArnon-Tirosh_Regulon_Boxplots_%s.pdf",grouping), device = "pdf",dpi = 300, width = 10, height = 6.68)
