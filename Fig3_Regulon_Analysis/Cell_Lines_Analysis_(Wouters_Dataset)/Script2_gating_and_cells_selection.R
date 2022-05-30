###################################################################################################
# (SCENIC) Gating signature enrichment AUCell Scores to determine phenotypes of 
# individual cells from Wouters et al scRNA-Seq data
# Natacha Comandante-Lou
# Nov, 2021
################################################################################################### 
rm(list=ls()) #clear all
cat("\014")  #clc
setwd("/~/Code/Fig3")

library(ggsci)
library(dplyr)
library(rstatix)
library(stringr)
library(ggpubr)
library(ggbeeswarm)
# Read Data --------------------------------------------------------------
grouping = "M-MT_MT-T-TN_TN-N-NU_NU-U"
X = read.csv(sprintf("Wouters_SCENIC_Merged_Diff_Reg_AUCell_DF_%s.csv",grouping)) #Supplemental Data 
colnames(X) = str_replace(colnames(X),'.x','')%>%str_replace('.y','-1')%>%str_replace('.3','-1')
colnames(X) = str_replace(colnames(X),'\\.\\.',' (')%>%str_replace('\\.',')')
X$cell_line = as.factor(X$cell_line)

regulon_list = c("FOS_motif_regulon","FOSB_motif_regulon","FOSL1_motif_regulon","FOSL2_motif_regulon",
                 "ATF1_motif_regulon","ATF2_motif_regulon","ATF4_motif_regulon",
                 "ATF5_motif_regulon","ATF6_motif_regulon","ATF6B_motif_regulon","ATF7_motif_regulon",
                 "JUN_motif_regulon","JUNB_motif_regulon","JUND_motif_regulon")
x = select(X,c("cell_id","cell_line","M_MT (187g)","MT_T_TN (179g)","TN_N_NU (197g)","NU_U (224g)",regulon_list))

#thresholds for each differentiation state
#1 - M_MT (187g)
#2 - MT_T_TN (179g)
#3 - TN_N_NU (197g)
#4 - NU_U (224g)

#melanocytic thresholds
thres.m.1 = 0.1
thres.m.3 = 0.02
thres.m.4 = 0.02

#transitory thresholds
thres.t.2 = 0.06
thres.t.4 = 0.02

#neural crest-like thresholds
thres.n.1 = 0.02
thres.n.3 = 0.05

#undifferentiated threshold
thres.u.1 = 0.02
thres.u.2 = 0.02
thres.u.4 = 0.06
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
count.diff.cellline = group_by(x,Diff_class, cell_line)%>%tally()


pal1 = pal_d3(palette= "category20")(12)
pal2 = pal_rickandmorty(palette = c("schwifty"), alpha = 0.9)(12)
mypal = c(pal1,pal2)

x$Diff_class = factor(x$Diff_class,levels = rev(c("Melanocytic","Transitory","Neural crest-like","Undifferentiated","NA")))


# Box plots (t-test) --------------------

.df = filter(x,!Diff_class%in%"NA")%>%group_by(Diff_class)
.df$Diff_class = factor(.df$Diff_class, levels = c("Melanocytic","Transitory","Neural crest-like","Undifferentiated"))

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
    geom_quasirandom(groupOnX=TRUE,size = 1,aes(color = cell_line))+
    geom_boxplot(alpha = 0.3,outlier.shape = NA)+
    scale_color_manual(values = mypal[sort(as.numeric(unique(.df[["cell_line"]])),decreasing = TRUE)])+
    guides(color = guide_legend(override.aes = list(size=10),title = "cell_line"))+
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

figs = lapply(c("FOSL2_motif_regulon","JUN_motif_regulon","FOSL1_motif_regulon","FOS_motif_regulon"),reg_boxplots_func,.df)

ggarrange(plotlist = figs, nrow = 1, ncol = 5, common.legend = TRUE)

ggsave(sprintf("Wouters_Regulon_Boxplots_%s.pdf",grouping), device = "pdf",dpi = 300, width = 12.55, height = 6.68)


