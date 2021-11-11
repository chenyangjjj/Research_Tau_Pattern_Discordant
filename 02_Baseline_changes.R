## This is for comparisions Comparisons of cortical tau deposition among different CSF/PET groups. 

graphics.off()
rm(list=ls())

library(stringr)
library(dplyr)
library(rprojroot)
library(ggseg)
library(ggpubr)
library(sjPlot) # For tab_model
library(doBy)
########################################################################################################## 
## Cutoffs of AD biomarkers based on ROC negative CN Vs. positive AD + MCI
Cutoff_CSF_ABETA42 = 978
Cutoff_CSF_PTAU = 23
Cutoff_CSF_TAU = 234
Cutoff_CSF_ABETA42.ABETA40 = 0.054
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
##################################################################################################################
## Still cutoffs form LONI and GMM method for ABETA42/40 ratio
Cutoff_CSF_MASS_ABETA42 = 1079
Cutoff_CSF_MASS_ABETA42.ABETA40 = 0.138
Cutoff_CSF_MASS_ABETA42.ABETA38 = 0.627
##################################################################################################################
Cutoff_SUVR_CER = 1.11
Cutoff_SUVR_Big_Ref = 0.82
Cutoff_SUVR_FBB_CER = 1.08
Cutoff_Centiloid = 24.4
Col_CU_CI = c("royalblue","darkorange")
Col_2 = c("#377EB8", "#E41A1C")
Col_3 = c("#4DAF4A", "#377EB8", "#E41A1C")
Col_4 = c("#4DAF4A", "#377EB8","coral1", "#E41A1C")
Col_5 = c("#4DAF4A", "#377EB8","lightpink1","coral1", "#E41A1C")
Col_8 = c("limegreen","royalblue2","steelblue1","blue3",
          "lightpink1","purple1","salmon1","red3") 
Data_ADNI_CSF_PET_Abeta_Tau = 
  read.csv("D:/Project/guo_20210206/CSF_PET_tau/Data/Data/Data_ADNI_AV1451_CSF_PET_Abeta_Tau_06_30_21.csv", 
           header=TRUE, sep=",",stringsAsFactors=FALSE)
Data_ADNI_CSF_PET_Abeta_Tau = 
  Data_ADNI_CSF_PET_Abeta_Tau[,2:ncol(Data_ADNI_CSF_PET_Abeta_Tau)]
Data_ADNI_CSF_PET_Abeta_Tau$Diag_AV1451_v1 = 
  factor(Data_ADNI_CSF_PET_Abeta_Tau$Diag_AV1451_v1 ,levels = c("CU","MCI","AD"))

# set new cut-off; 
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta=
  factor(Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta,
         levels = c("CSF-/PET-","CSF+/PET-","CSF-/PET+","CSF+/PET+"))

Data_ADNI_CSF_PET_Abeta_Tau_v1 = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau, AV1451_Scan_ID==1)
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1 = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau_v1, 
         (Diag_AV1451_v1 != "AD"))

Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1, 
         !((Stage_CSF_PET_Abeta == "CSF-/PET-")&
             ((Positivity_CSF_PTAU_ABETA40_Closest_AV1451_v1 == "CSF_PTAU_ABETA40_P")|
                (Positivity_AV1451_ENTORHINAL_Non_PVC_Inf_CER_GM_v1 == "FTP_P")|
                (Positivity_AV1451_Meta_TEM_ROI_Non_PVC_Inf_CER_GM_v1 == "FTP_P"))))
nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
################################################################################################################## 去掉两个netorhinal tau PET 异常被试

# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A =
#   subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
#          !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 1.6)&(Stage_CSF_PET_Abeta == "CSF+/PET-")))
# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 =
#   subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A,
#          !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 2)&(Stage_CSF_PET_Abeta == "CSF+/PET+")))

Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[is.na(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$
                                                         APOE_status),c("APOE_status")] = "APOE4_Non_carrier"
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage)

Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref = 
subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, Stage_CSF_PET_Abeta == "CSF-/PET-")

## CSF_ABETA42.ABETA40_Closest_AV1451_v1
################################################################################################################## 
# CSF_ABETA42.ABETA40_Closest_AV1451_v1
GLM_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta_Vs_Ref= 
  glm(CSF_ABETA42.ABETA40_Closest_AV1451_v1~Stage_CSF_PET_Abeta+
        Gender + Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
tab_model(GLM_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta_Vs_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")

Median_CSF_ABETA42_ABETA40_Ref = 
  median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1)

Plot_Comparison_BL_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         aes(x = Stage_CSF_PET_Abeta, y = CSF_ABETA42.ABETA40_Closest_AV1451_v1)) +
  geom_point(aes(color=Stage_CSF_PET_Abeta,
                 fill=Stage_CSF_PET_Abeta), size=1, shape=19,
             position=position_jitter(width=0.25, height=0),alpha = 0.5) +
  geom_boxplot(outlier.colour=NA, fill=NA, aes(colour=Stage_CSF_PET_Abeta), fatten=0.7) +
  scale_color_manual(values = Col_4)+
  annotate("text", x =  1, y = 0.125,
           label = 'Ref',
           size = 2.5, color = Col_4[1]) +
  annotate("text", x =  2, y = 0.06,
           label = expression(bold(italic(p)*' < 0.001')),
           size = 2.5, color = Col_4[2]) +
  annotate("text", x =  3, y = 0.1125,
           label = expression(bold(italic(p)*' = 0.0003')),
           size = 2.5,color = Col_4[3]) +
  annotate("text", x =  4, y = 0.06,
           label = expression(bold(italic(p)*' < 0.001')),
           size = 2.5,color = Col_4[4]) +
  scale_y_continuous(breaks = seq(0.025, 0.125,0.025),limits = c(0.015, 0.125)) +
  xlab("") +
  ylab(expression('CSF A'*beta[42]*'/A'*beta[40])) +
  geom_hline(aes(yintercept=Median_CSF_ABETA42_ABETA40_Ref),
             color=Col_4[1],linetype="dotted", size=0.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(title=element_text(size=9),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=8),axis.title=element_text(size=9,face="plain"),
        legend.position = "none", legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
# ggsave("Plot_Comparison_BL_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta.jpg",
#        Plot_Comparison_BL_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta,
#        units="in", width=5, height=5, dpi=600)
# ggsave("Plot_Comparison_BL_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta.pdf",
#        Plot_Comparison_BL_CSF_ABETA42_ABETA40_Stage_CSF_PET_Abeta,
#        units="in", width=5, height=5)


## PET_ABETA_CL_Closest_AV1451_v1
####################################################################################################
# PET_ABETA_CL_Closest_AV1451_v1
summaryBy(PET_ABETA_CL_Closest_AV1451_v1 ~ Stage_CSF_PET_Abeta, 
          data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
          FUN = list(median,IQR))
compare_means(PET_ABETA_CL_Closest_AV1451_v1 ~ Stage_CSF_PET_Abeta,
              data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
              method = "wilcox.test", paired = FALSE,p.adjust.method ="BH")
Median_PET_ABETA_CL_Ref = 
  median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref$PET_ABETA_CL_Closest_AV1451_v1)
GLM_PET_ABETA_CL_Stage_CSF_PET_Abeta_Vs_Ref= 
  glm(PET_ABETA_CL_Closest_AV1451_v1~Stage_CSF_PET_Abeta+
        Gender + Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
tab_model(GLM_PET_ABETA_CL_Stage_CSF_PET_Abeta_Vs_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
Plot_Comparison_BL_PET_ABETA_CL_Stage_CSF_PET_Abeta = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         aes(x = Stage_CSF_PET_Abeta, y = PET_ABETA_CL_Closest_AV1451_v1)) +
  geom_point(aes(color=Stage_CSF_PET_Abeta,
                 fill=Stage_CSF_PET_Abeta), size=1, shape=19,
             position=position_jitter(width=0.25, height=0),alpha = 0.5) +
  geom_boxplot(outlier.colour=NA, fill=NA, aes(colour=Stage_CSF_PET_Abeta), fatten=0.7) +
  scale_color_manual(values = Col_4)+
  annotate("text", x =  1, y = 45,
           label = 'Ref',
           size = 2.5, color = Col_4[1]) +
  annotate("text", x =  2, y = 45,
           label = expression(italic(p)*' = 0.23'),
           size = 2.5, color = Col_4[2]) +
  annotate("text", x =  3, y = 70,
           label = expression(bold(italic(p)*' < 0.001')),
           size = 2.5,color = Col_4[3]) +
  annotate("text", x =  4, y = 160,
           label = expression(bold(italic(p)*' < 0.001')),
           size = 2.5,color = Col_4[4]) +
  scale_y_continuous(breaks = seq(-50, 250,50),limits = c(-50, 250)) +
  xlab("") +
  ylab(expression('A'*beta*' PET (Centiloids)')) +
  geom_hline(aes(yintercept=Median_PET_ABETA_CL_Ref),
             color=Col_4[1],linetype="dotted", size=0.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(title=element_text(size=9),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=8),axis.title=element_text(size=11,face="plain"),
        legend.position = "none", legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
# ggsave("Plot_Comparison_BL_PET_ABETA_CL_Stage_CSF_PET_Abeta.jpg",
#        Plot_Comparison_BL_PET_ABETA_CL_Stage_CSF_PET_Abeta,
#        units="in", width=5, height=5, dpi=600)
# ggsave("Plot_Comparison_BL_PET_ABETA_CL_Stage_CSF_PET_Abeta.pdf",
#        Plot_Comparison_BL_PET_ABETA_CL_Stage_CSF_PET_Abeta,
#        units="in", width=5, height=5)


## CSF_PTAU.ABETA40_Closest_AV1451_v1
####################################################################################################
# CSF_PTAU.ABETA40_Closest_AV1451_v1
summaryBy(CSF_PTAU.ABETA40_Closest_AV1451_v1 ~ Stage_CSF_PET_Abeta, 
          data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
          FUN = list(median,IQR))
compare_means(CSF_PTAU.ABETA40_Closest_AV1451_v1 ~ Stage_CSF_PET_Abeta,
              data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
              method = "wilcox.test", paired = FALSE,p.adjust.method ="BH")
Median_CSF_PTAU_ABETA40_Ref = 
  median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1)
GLM_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta_Vs_Ref= 
  glm(CSF_PTAU.ABETA40_Closest_AV1451_v1~Stage_CSF_PET_Abeta+
        Gender + Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
tab_model(GLM_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta_Vs_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
Plot_Comparison_BL_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         aes(x = Stage_CSF_PET_Abeta, y = CSF_PTAU.ABETA40_Closest_AV1451_v1)) +
  geom_point(aes(color=Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1,
                 fill=Stage_CSF_PET_Abeta), size=1, shape=19,
             position=position_jitter(width=0.25, height=0),alpha = 0.5) +
  geom_boxplot(outlier.colour=NA, fill=NA, aes(colour=Stage_CSF_PET_Abeta), fatten=0.7) +
  scale_color_manual(values = Col_4)+
  scale_shape_manual(values = c(19,1)) +
  annotate("text", x =  1, y = 0.0016,
           label = 'Ref',
           size = 2.5, color = Col_4[1]) +
  annotate("text", x =  2, y = 0.002,
           label = expression(bold(italic(p)*' = 0.002')),
           size = 2.5, color = Col_4[2]) +
  annotate("text", x =  3, y = 0.0016,
           label = expression(italic(p)*' = 0.38'),
           size = 2.5,color = Col_4[3]) +
  annotate("text", x =  4, y = 0.0037,
           label = expression(bold(italic(p)*' < 0.001')),
           size = 2.5,color = Col_4[4]) +
  scale_y_continuous(breaks = seq(0.001, 0.004,0.001),limits = c(0.0005, 0.0046)) +
  xlab("") +
  ylab(expression('CSF p-Tau/A'*beta[40])) +
  geom_hline(aes(yintercept=Median_CSF_PTAU_ABETA40_Ref),
             color=Col_4[1],linetype="dotted", size=0.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(title=element_text(size=9),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=8),axis.title=element_text(size=11,face="plain"),
        legend.position = "right", legend.direction = "vertical")+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
# ggsave("Plot_Comparison_BL_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta.jpg",
#        Plot_Comparison_BL_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5, dpi=600)
# ggsave("Plot_Comparison_BL_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta.pdf",
#        Plot_Comparison_BL_CSF_PTAU_ABETA40_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5)


## ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1

# ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1
####################################################################################################
summaryBy(ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ Stage_CSF_PET_Abeta, 
          data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
          FUN = list(median,IQR))
Median_ENTORHINAL_FTP_Ref = 
  median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
GLM_ENTORHINAL_FTP_Stage_CSF_PET_Abeta_Vs_Ref= 
  glm(ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~Stage_CSF_PET_Abeta+
        Gender + Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
tab_model(GLM_ENTORHINAL_FTP_Stage_CSF_PET_Abeta_Vs_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
Plot_Comparison_BL_ENTORHINAL_FTP_Stage_CSF_PET_Abeta = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         aes(x = Stage_CSF_PET_Abeta, y = ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_point(aes(color=Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1,
                 fill=Stage_CSF_PET_Abeta), size=1.5,
             position=position_jitter(width=0.25, height=0),alpha = 0.7) +
  geom_boxplot(outlier.colour=NA, fill=NA, aes(colour=Stage_CSF_PET_Abeta), fatten=0.7) +
  scale_color_manual(values = Col_4)+
  scale_shape_manual(values = c(19,1)) +
  annotate("text", x =  1, y = 1.4,
           label = 'Ref',
           size = 2.5, color = Col_4[1]) +
  annotate("text", x =  2, y = 1.4,
           label = expression(italic(P)*' = 0.012'),
           size = 2.5, color = Col_4[2]) +
  annotate("text", x =  3, y = 1.4,
           label = expression(italic(P)*' = 0.19'),
           size = 2.5,color = Col_4[3]) +
  annotate("text", x =  4, y = 2.0,
           label = expression(bold(italic(P)*' < 0.001')),
           size = 2.5,color = Col_4[4]) +
  scale_y_continuous(breaks = seq(0.8, 2.4,0.4),limits = c(0.8, 2.4)) +
  xlab("") +
  ylab("Entorhinal FTP SUVR") +
  geom_hline(aes(yintercept=Median_ENTORHINAL_FTP_Ref),
             color=Col_4[1],linetype="dotted", size=0.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(title=element_text(size=9),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=8),axis.title=element_text(size=11,face="plain"),
        legend.position = "none", legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

stage4 = subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,Stage_CSF_PET_Abeta == "CSF+/PET-")
hist(stage4$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)

# ggsave("Plot_Comparison_BL_ENTORHINAL_FTP_Stage_CSF_PET_Abeta.jpg",
#        Plot_Comparison_BL_ENTORHINAL_FTP_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5, dpi=600)
# ggsave("Plot_Comparison_BL_ENTORHINAL_FTP_Stage_CSF_PET_Abeta.pdf",
#        Plot_Comparison_BL_ENTORHINAL_FTP_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5)




## BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1

####################################################################################################
# BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1
####################################################################################################
summaryBy(BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ Stage_CSF_PET_Abeta, 
          data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
          FUN = list(median,IQR))
compare_means(BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ Stage_CSF_PET_Abeta,
              data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
              method = "wilcox.test", paired = FALSE,p.adjust.method ="BH")
Median_BRAAK34_FTP_Ref = 
  median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
GLM_BRAAK34_FTP_Stage_CSF_PET_Abeta_Vs_Ref= 
  glm(BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~Stage_CSF_PET_Abeta+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
tab_model(GLM_BRAAK34_FTP_Stage_CSF_PET_Abeta_Vs_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
Plot_Comparison_BL_BRAAK34_FTP_Stage_CSF_PET_Abeta = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         aes(x = Stage_CSF_PET_Abeta, y = BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_point(aes(color=Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1,
                 fill=Stage_CSF_PET_Abeta), size=1.5,
             position=position_jitter(width=0.25, height=0),alpha = 0.7) +
  geom_boxplot(outlier.colour=NA, fill=NA, aes(colour=Stage_CSF_PET_Abeta), fatten=0.7) +
  scale_color_manual(values = Col_4)+
  scale_shape_manual(values = c(19,1)) +
  annotate("text", x =  1, y = 1.4,
           label = 'Ref',
           size = 2.5, color = Col_4[1]) +
  annotate("text", x =  2, y = 1.4,
           label = expression(italic(P)*' = 0.21'),
           size = 2.5, color = Col_4[2]) +
  annotate("text", x =  3, y = 1.4,
           label = expression(italic(P)*' = 0.07'),
           size = 2.5,color = Col_4[3]) +
  annotate("text", x =  4, y = 2.2,
           label = expression(bold(italic(P)*' < 0.001')),
           size = 2.5,color = Col_4[4]) +
  scale_y_continuous(breaks = seq(0.8,2.4,0.4),limits = c(0.8,2.4)) +
  xlab("") +
  labs(y=expression(paste(Braak[III/IV],' FTP SUVR')))+
  geom_hline(aes(yintercept = Median_BRAAK34_FTP_Ref),
             color=Col_4[1],linetype="dotted", size=0.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(title=element_text(size=9),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=8),axis.title=element_text(size=11,face="plain"),
        legend.position = "none", legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# ggsave("Plot_Comparison_BL_BRAAK34_FTP_Stage_CSF_PET_Abeta.jpg",
#        Plot_Comparison_BL_BRAAK34_FTP_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5, dpi=600)
# ggsave("Plot_Comparison_BL_BRAAK34_FTP_Stage_CSF_PET_Abeta.pdf",
#        Plot_Comparison_BL_BRAAK34_FTP_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5)

stage4 = subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,Stage_CSF_PET_Abeta == "CSF+/PET+")
stage4$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1

# BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1
####################################################################################################
summaryBy(BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ Stage_CSF_PET_Abeta, 
          data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
          FUN = list(median,IQR))
compare_means(BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ Stage_CSF_PET_Abeta,
              data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
              method = "wilcox.test", paired = FALSE,p.adjust.method ="BH")
Median_BRAAK56_FTP_Ref = 
  median(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
GLM_BRAAK56_FTP_Stage_CSF_PET_Abeta_Vs_Ref= 
  glm(BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~Stage_CSF_PET_Abeta+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
tab_model(GLM_BRAAK56_FTP_Stage_CSF_PET_Abeta_Vs_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         aes(x = Stage_CSF_PET_Abeta, y = BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_point(aes(color=Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1,
                 fill=Stage_CSF_PET_Abeta), size=1.5,
             position=position_jitter(width=0.25, height=0),alpha = 0.7) +
  geom_boxplot(outlier.colour=NA, fill=NA, aes(colour=Stage_CSF_PET_Abeta), fatten=1) +
  scale_color_manual(values = Col_4)+
  scale_shape_manual(values = c(19,1)) +
  annotate("text", x =  1, y = 1.4,
           label = 'Ref',
           size = 2.5, color = Col_4[1]) +
  annotate("text", x =  2, y = 1.4,
           label = expression(italic(P)*' = 0.39'),
           size = 2.5, color = Col_4[2]) +
  annotate("text", x =  3, y = 1.4,
           label = expression(bold(italic(P)*' = 0.004')),
           size = 2.5,color = Col_4[3]) +
  annotate("text", x =  4, y = 1.9,
           label = expression(bold(italic(P)*' < 0.001')),
           size = 2.5,color = Col_4[4]) +
  scale_y_continuous(breaks = seq(0.8,2,0.4),limits = c(0.8,2)) +
  xlab("") +
  labs(y=expression(paste(Braak[V/VI],' FTP SUVR')))+
  geom_hline(aes(yintercept = Median_BRAAK56_FTP_Ref),
             color=Col_4[1],linetype="dotted", size=0.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  theme(title=element_text(size=9),plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=8),axis.title=element_text(size=11,face="plain"),
        legend.position = "none", legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta
# ggsave("Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta.jpg",
#        Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5, dpi=600)
# ggsave("Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta.pdf",
#        Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta, 
#        units="in", width=5, height=5)

Plot_Comparison_BL_ABETA_Tau_Stages_CSF_PET_Abeta = 
  ggarrange(
    Plot_Comparison_BL_ENTORHINAL_FTP_Stage_CSF_PET_Abeta,
    Plot_Comparison_BL_BRAAK34_FTP_Stage_CSF_PET_Abeta,
    Plot_Comparison_BL_BRAAK56_FTP_Stage_CSF_PET_Abeta,
    nrow=1,ncol = 3, labels = c("D", "E","F"),
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE,legend = "bottom")


ggsave("D:/Project/guo_20210206/CSF_PET_tau/Output/main_results/with_otuliers/1_with_outliers/F2_new.tiff",
       Plot_Comparison_BL_ABETA_Tau_Stages_CSF_PET_Abeta,
       units="in", width=10.2, height=4, dpi=300)






