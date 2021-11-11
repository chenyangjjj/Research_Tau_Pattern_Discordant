## This code is to find different corresponding tau pattern of different CSF/PET groups categorized 
## by CSF Abeta and Aβ PET
## main part

graphics.off()
rm(list=ls())
library(stringr)
library(dplyr)
library(rprojroot)
library(ggseg)
library(ggpubr)
library(tidyverse)

########################################################################################################## 
## Cutoffs of AD biomarkers based on ROC negative CN Vs. positive AD + MCI
########################################################################################################## 
Cutoff_CSF_ABETA42 = 978
Cutoff_CSF_PTAU = 23
Cutoff_CSF_TAU = 234
Cutoff_CSF_ABETA42.ABETA40 = 0.054  # 0.051 for Gaussian-mixture model; 0.054 for ROC analysis youden index
Cutoff_CSF_PTAU.ABETA40 = 0.0012
Cutoff_aHCV = 6787
Cutoff_Meta_ROI_thickness = 2.6
Cutoff_FDG_Ref_PON = 1.21
Cutoff_FTP_Inf_CER_GM_PVC_Meta_TEM_ROI_SUVR = 1.55
Cutoff_FTP_Inf_CER_GM_Non_PVC_Meta_TEM_ROI_SUVR = 1.25
Cutoff_FTP_Inf_CER_GM_PVC_ENTORHINAL_SUVR = 1.90
Cutoff_FTP_Inf_CER_GM_Non_PVC_ENTORHINAL_SUVR = 1.21
Cutoff_FTP_WM_PVC_Meta_TEM_ROI_SUVR = 2.05
Cutoff_FTP_WM_Non_PVC_Meta_TEM_ROI_SUVR = 1.07
Cutoff_FTP_WM_PVC_ENTORHINAL_SUVR = 2.23
Cutoff_FTP_WM_Non_PVC_ENTORHINAL_SUVR = 1.06
## Still cutoffs form LONI and GMM method for ABETA42/40 ratio
Cutoff_CSF_MASS_ABETA42 = 1079
Cutoff_CSF_MASS_ABETA42.ABETA40 = 0.138
Cutoff_CSF_MASS_ABETA42.ABETA38 = 0.627
Cutoff_SUVR_CER = 1.11
Cutoff_SUVR_Big_Ref = 0.82
Cutoff_SUVR_FBB_CER = 1.08
Cutoff_Centiloid = 24.4  # this is the new Aβ PET cut-off
Data_ADNI_CSF_PET_Abeta_Tau = 
  read.csv("D:/Project/guo_20210206/CSF_PET_tau/Data/Data/Data_ADNI_AV1451_CSF_PET_Abeta_Tau_06_30_21.csv", 
           header=TRUE, sep=",",stringsAsFactors=FALSE)
Data_ADNI_CSF_PET_Abeta_Tau = 
  Data_ADNI_CSF_PET_Abeta_Tau[,2:ncol(Data_ADNI_CSF_PET_Abeta_Tau)]
Data_ADNI_CSF_PET_Abeta_Tau$Diag_AV1451_v1 = 
  factor(Data_ADNI_CSF_PET_Abeta_Tau$Diag_AV1451_v1 ,levels = c("CU","MCI","AD"))

########################################################################################################## 
# define stages  CSF-/PET-  CSF+/PET-  CSF-/PET+  CSF+/PET+
########################################################################################################## 
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta=
  factor(Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta,
         levels = c("CSF-/PET-","CSF+/PET-","CSF-/PET+","CSF+/PET+"))


########################################################################################################## 
# remove AD participants and FTP SUVR outliers
########################################################################################################## 
Data_ADNI_CSF_PET_Abeta_Tau_v1 = subset(Data_ADNI_CSF_PET_Abeta_Tau, AV1451_Scan_ID==1)
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1 = subset(Data_ADNI_CSF_PET_Abeta_Tau_v1, (Diag_AV1451_v1 != "AD"))
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1$Stage_CSF_PET_Abeta)
nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1)
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1, 
         !((Stage_CSF_PET_Abeta == "CSF-/PET-")&
             ((Positivity_CSF_PTAU_ABETA40_Closest_AV1451_v1 == "CSF_PTAU_ABETA40_P")|
                (Positivity_AV1451_ENTORHINAL_Non_PVC_Inf_CER_GM_v1 == "FTP_P")|
                (Positivity_AV1451_Meta_TEM_ROI_Non_PVC_Inf_CER_GM_v1 == "FTP_P"))))
# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A =
#   subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
#          !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 1.6)&(Stage_CSF_PET_Abeta == "CSF+/PET-")))
# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 =
#   subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A,
#          !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 2)&(Stage_CSF_PET_Abeta == "CSF+/PET+")))

Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[is.na(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$
                                                         APOE_status),c("APOE_status")] = "APOE4_Non_carrier"
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta)



# 先画qqplot,目视检查，再结合显著性检验，若p < 0.05, 则非正态
# shapiro test检验的是数据的分布与正态分布之间有无显著差异
shapiro.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
ggqqplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
       aes(x=ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) + 
  geom_histogram(color="gray50",bins=25)+
  xlab("") +
  ylab("Values") +
  guides(colour = guide_legend(nrow =1)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=9),
        plot.title = element_text(size=9,face="plain",hjust = 0.5),
        axis.text=element_text(size=9),axis.title=element_text(size=9,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


Data_Label_BL_AV1451_Inf_CER_GM_Non_PVC = 
  Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1[, grepl("CTX_",
                                                names( Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1))&
                                          grepl("_Inf_CER_GM_Non_PVC_v1", 
                                                names( Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1))]
Region_DTK_BL_AV1451_Inf_CER_GM_Non_PVC = colnames(Data_Label_BL_AV1451_Inf_CER_GM_Non_PVC)
DTK_Labes_For_Plot = 
  c("bankssts","caudal anterior cingulate","caudal middle frontal",
    "cuneus","entorhinal","frontal pole","fusiform","inferior parietal",
    "inferior temporal","insula","isthmus cingulate","lateral occipital",
    "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal",
    "paracentral","parahippocampal","pars opercularis","pars orbitalis",
    "pars triangularis","pericalcarine","postcentral","posterior cingulate",
    "precentral","precuneus","rostral anterior cingulate","rostral middle frontal",
    "superior frontal","superior parietal","superior temporal","supramarginal",
    "temporal pole","transverse temporal")
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_L = as.data.frame(DTK_Labes_For_Plot)
colnames(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_L)[1] = "region"
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_L$hemi="left"
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_R = as.data.frame(DTK_Labes_For_Plot)
colnames(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_R)[1] = "region"
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_R$hemi="right"
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta = 
  rbind(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_L,Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta_R)
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label = Region_DTK_BL_AV1451_Inf_CER_GM_Non_PVC
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage1=NA
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage2=NA
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage3=NA
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage4=NA
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage2_Vs_Stage1=NA
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage3_Vs_Stage1=NA
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage4_Vs_Stage1=NA
attach(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta)

for (i in 1:68)
{
  #i=1
  Region_AV1451_v1 = Region_DTK_BL_AV1451_Inf_CER_GM_Non_PVC[i]
  # Region_AV1451_Slope = Region_DTK_Slope_AV1451_WM_Non_PVC[i]
  ##########################################################################################################   
  ## Mean AV1451 SUVR in different ROIs for different A/T/N profiles
  ########################################################################################################## 
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[,Region_AV1451_v1][Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta
                                                                               =="CSF-/PET-"])
  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage2[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[,Region_AV1451_v1][Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta
                                                                               =="CSF+/PET-" ])
  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage3[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[,Region_AV1451_v1][Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta
                                                                               =="CSF-/PET+"])
  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$Baseline_AV1451_Stage4[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[,Region_AV1451_v1][Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta
                                                                               =="CSF+/PET+"])
  
  
  ##########################################################################################################   
  ## compare AV1451 SUVR with reference CSF-/PET- in different ROIs for different CSF/PET groups
  ##########################################################################################################  
  # funtion_glm = paste( Region_AV1451_v1,"~Stage_CSF_PET_Abeta+Diag_AV1451_v1+Gender+Age_AV1451_v1+APOE_status")
  
  funtion_glm = paste(Region_AV1451_v1,"~Stage_CSF_PET_Abeta+Gender+Age_AV1451_v1")
  GLM_Regional_AV1451_Stage_CSF_PET_Abeta = 
    glm(funtion_glm,family = gaussian,  
        data =Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
  
  Sta_GLM_Regional_AV1451_Stage_CSF_PET_Abeta= summary(GLM_Regional_AV1451_Stage_CSF_PET_Abeta)
  P_values_GLM = Sta_GLM_Regional_AV1451_Stage_CSF_PET_Abeta$coefficients[2:4,3:4]
  # View(P_values_GLM)
  Data_Sta_GLM = as.data.frame(P_values_GLM)
  colnames(Data_Sta_GLM)[1] = "T_values_GLM"
  colnames(Data_Sta_GLM)[2] = "P_values_GLM"
  Data_Sta_GLM$Stage_Label = str_sub(row.names(Data_Sta_GLM),-9)
  Data_Sta_GLM$Stage_Label_b = c("Stage2","Stage3","Stage4")
  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage2_Vs_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    Data_Sta_GLM$P_values_GLM[Data_Sta_GLM$Stage_Label_b=="Stage2"]
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage3_Vs_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    Data_Sta_GLM$P_values_GLM[Data_Sta_GLM$Stage_Label_b=="Stage3"]  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage4_Vs_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    Data_Sta_GLM$P_values_GLM[Data_Sta_GLM$Stage_Label_b=="Stage4"]
  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$t_AV1451_Stage2_Vs_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    Data_Sta_GLM$T_values_GLM[Data_Sta_GLM$Stage_Label_b=="Stage2"]
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$t_AV1451_Stage3_Vs_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    Data_Sta_GLM$T_values_GLM[Data_Sta_GLM$Stage_Label_b=="Stage3"]  
  Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$t_AV1451_Stage4_Vs_Stage1[Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$area_label==Region_AV1451_v1] = 
    Data_Sta_GLM$T_values_GLM[Data_Sta_GLM$Stage_Label_b=="Stage4"]
}

detach(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta)

Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage2_Vs_Stage1_BH  = 
  p.adjust(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage2_Vs_Stage1, method = "BH",
                 n = 68)

Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage3_Vs_Stage1_BH  = 
  p.adjust(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage3_Vs_Stage1, method = "BH",
           n = 68)

Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage4_Vs_Stage1_BH  = 
  p.adjust(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage4_Vs_Stage1, method = "BH",
           n = 68)

View(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta)
Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta$p_AV1451_Stage4_Vs_Stage1

#write.csv(x = sig_stage3,file = 'E:/Project01_CSF_pet/1_Results/Sig_stage3_vs_Ref.csv')

Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage1 = 
  ggseg(.data = Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta,atlas = "dk",
        mapping=aes(fill=Baseline_AV1451_Stage1), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","firebrick","goldenrod"),
                       breaks=seq(1,1.6,0.2),limits=c(1,1.6))+
  xlab("") +
  ylab("") +
  labs(title ="CSF-/PET-",fill="SUVR") +
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage1_Long_Legend = 
  ggseg(.data = Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta,atlas = "dk",
        mapping=aes(fill=Baseline_AV1451_Stage1), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","firebrick","goldenrod"),
                       breaks=seq(1,1.6,0.2),limits=c(1,1.6))+
  xlab("") +
  ylab("") +
  labs(title ="CSF-/PET-",fill="SUVR") +
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.y = unit(0.15, 'cm'),legend.key.width = unit(1.5, "cm"),legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


########################################################################################################## 
## CSF+/PET-   Baseline FTP SUVR
########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage2 = 
  ggseg(.data = Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta,atlas = "dk",
        mapping=aes(fill=Baseline_AV1451_Stage2), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","firebrick","goldenrod"),
                       breaks=seq(1,1.6,0.2),limits=c(1,1.6))+
  xlab("") +
  ylab("") +
  labs(title ="CSF+/PET-",fill="SUVR") +
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 



########################################################################################################## 
## CSF-/PET+   Baseline FTP SUVR
########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage3 = 
  ggseg(.data = Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta,atlas = "dk",
        mapping=aes(fill=Baseline_AV1451_Stage3), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","firebrick","goldenrod"),
                       breaks=seq(1,1.6,0.2),limits=c(1,1.6))+
  xlab("") +
  ylab("") +
  labs(title ="CSF-/PET+",fill="SUVR") +
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 





########################################################################################################## 
## Compare baseline FTP SUVR between CSF+/PET- and CSF-/PET-  without BH correction
########################################################################################################## 
Baseline_AV1451_Stage_CSF_PET_Abeta_Stage2_Vs_Stage1_Sig =
  subset(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta, (p_AV1451_Stage2_Vs_Stage1<0.05)&
           (t_AV1451_Stage2_Vs_Stage1>0))

########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage2_Vs_Stage1_Long_Legend = 
  ggseg(.data = Baseline_AV1451_Stage_CSF_PET_Abeta_Stage2_Vs_Stage1_Sig,atlas = "dk",
        mapping=aes(fill=p_AV1451_Stage2_Vs_Stage1), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("goldenrod","firebrick","light blue","dodgerblue4"),
                       breaks=seq(0,0.05,0.01),limits=c(0,0.05))+
  xlab("") +
  ylab("") +
  labs(title ="CSF+/PET- > CSF-/PET-",fill="p value") +
  # theme(legend.title=element_blank())+
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.x = unit(0.15, 'cm'),legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


########################################################################################################## 
## Significant higher tau than reference group with BH correction
########################################################################################################## 
Baseline_AV1451_Stage_CSF_PET_Abeta_Stage2_Vs_Stage1_Sig_BH =
  subset(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta, (p_AV1451_Stage2_Vs_Stage1_BH<0.05)&
           (t_AV1451_Stage2_Vs_Stage1>0))

########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage2_Vs_Stage1_BH_Long_Legend = 
  ggseg(.data = Baseline_AV1451_Stage_CSF_PET_Abeta_Stage2_Vs_Stage1_Sig_BH,atlas = "dk",
        mapping=aes(fill=p_AV1451_Stage2_Vs_Stage1_BH), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("goldenrod","firebrick","light blue","dodgerblue4"),
                       breaks=seq(0,0.05,0.01),limits=c(0,0.05))+
  xlab("") +
  ylab("") +
  labs(title ="CSF+/PET- > CSF-/PET-",fill="p value") +
  # theme(legend.title=element_blank())+
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.x = unit(0.15, 'cm'),legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


########################################################################################################## 
## Compare baseline FTP SUVR between CSF-/PET+ and CSF-/PET-   without BH correction
########################################################################################################## 
Baseline_AV1451_Stage_CSF_PET_Abeta_Stage3_Vs_Stage1_Sig =
  subset(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta, (p_AV1451_Stage3_Vs_Stage1<0.05)&
           (t_AV1451_Stage3_Vs_Stage1>0))
########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage3_Vs_Stage1 = 
  ggseg(.data = Baseline_AV1451_Stage_CSF_PET_Abeta_Stage3_Vs_Stage1_Sig,atlas = "dk",
        mapping=aes(fill=p_AV1451_Stage3_Vs_Stage1), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("goldenrod","firebrick","light blue","dodgerblue4"),
                       breaks=seq(0,0.05,0.01),limits=c(0,0.05))+
  xlab("") +
  ylab("") +
  labs(title ="CSF-/PET+ > CSF-/PET-",fill="p value") +
  # theme(legend.title=element_blank())+
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.x = unit(0.15, 'cm'),legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 



########################################################################################################## 
## Significant higher tau than reference group with BH correction
########################################################################################################## 
Baseline_AV1451_Stage_CSF_PET_Abeta_Stage3_Vs_Stage1_Sig_BH =
  subset(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta, (p_AV1451_Stage3_Vs_Stage1_BH<0.05)&
           (t_AV1451_Stage3_Vs_Stage1>0))
########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage3_Vs_Stage1_BH = 
  ggseg(.data = Baseline_AV1451_Stage_CSF_PET_Abeta_Stage3_Vs_Stage1_Sig_BH,atlas = "dk",
        mapping=aes(fill=p_AV1451_Stage3_Vs_Stage1_BH), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("goldenrod","firebrick","light blue","dodgerblue4"),
                       breaks=seq(0,0.05,0.01),limits=c(0,0.05))+
  xlab("") +
  ylab("") +
  labs(title ="CSF-/PET+ > CSF-/PET-",fill="p value") +
  # theme(legend.title=element_blank())+
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.x = unit(0.15, 'cm'),legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 



########################################################################################################## 
## Compare baseline FTP SUVR between CSF+/PET+ and CSF-/PET-   without BH correction
########################################################################################################## 
Baseline_AV1451_Stage_CSF_PET_Abeta_Stage4_Vs_Stage1_Sig =
  subset(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta, (p_AV1451_Stage4_Vs_Stage1<0.05)&
           (t_AV1451_Stage4_Vs_Stage1>0))
########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage4_Vs_Stage1 = 
  ggseg(.data = Baseline_AV1451_Stage_CSF_PET_Abeta_Stage4_Vs_Stage1_Sig,atlas = "dk",
        mapping=aes(fill=p_AV1451_Stage4_Vs_Stage1), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("goldenrod","firebrick","light blue","dodgerblue4"),
                       breaks=seq(0,0.05,0.01),limits=c(0,0.05))+
  xlab("") +
  ylab("") +
  labs(title ="CSF+/PET+ > CSF-/PET-",fill="p value") +
  # theme(legend.title=element_blank())+
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.x = unit(0.15, 'cm'),legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 



########################################################################################################## 
## Significant higher tau than reference group with BH correction
########################################################################################################## 
Baseline_AV1451_Stage_CSF_PET_Abeta_Stage4_Vs_Stage1_Sig_BH =
  subset(Median_Regional_AV1451_DTK_Stage_CSF_PET_Abeta, (p_AV1451_Stage4_Vs_Stage1_BH<0.05)&
           (t_AV1451_Stage4_Vs_Stage1>0))
########################################################################################################## 
Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage4_Vs_Stage1_BH = 
  ggseg(.data = Baseline_AV1451_Stage_CSF_PET_Abeta_Stage4_Vs_Stage1_Sig_BH,atlas = "dk",
        mapping=aes(fill=p_AV1451_Stage4_Vs_Stage1_BH), position="dispersed",
        colour="black",size=.1,show.legend = T) +
  scale_fill_gradientn(colours=c("goldenrod","firebrick","light blue","dodgerblue4"),
                       breaks=seq(0,0.05,0.01),limits=c(0,0.05))+
  xlab("") +
  ylab("") +
  labs(title ="CSF+/PET+ > CSF-/PET-",fill="p value") +
  # theme(legend.title=element_blank())+
  theme(title=element_text(size=12),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="plain"),
        legend.text=element_text(size=12),legend.position = "bottom", 
        legend.spacing.x = unit(0.15, 'cm'),legend.key.width = unit(2, "cm"),
        legend.direction = "horizontal")+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


############################################################################################
############################################################################################
mainDir = find_rstudio_root_file() 
setwd(file.path(mainDir,"/Output"))
dir.create(file.path("Spatial_pattern_Tau"), showWarnings = FALSE)
setwd(file.path("Spatial_pattern_Tau"))
############################################################################################
############################################################################################
############################################################################################

Plot_Comparison_Baseline_ROI_FTP_SUVR_Stages_CSF_PET_Abeta = 
  ggarrange(
    Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage2_Vs_Stage1_Long_Legend,
    Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage3_Vs_Stage1_BH,
    Plot_Comparison_Baseline_ROI_FTP_SUVR_Stage4_Vs_Stage1_BH,
    nrow=3,ncol = 1, labels = c("D", "E","F"),
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend="bottom")


ggsave("D:/Project/guo_20210206/CSF_PET_tau/Output/main_results/with_otuliers/1_with_outliers/Figure1_tau_pattern.tiff",
       Plot_Comparison_Baseline_ROI_FTP_SUVR_Stages_CSF_PET_Abeta,
       units="in", width=6.5, height=5, dpi=600)







