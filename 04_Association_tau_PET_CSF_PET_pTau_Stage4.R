## This code is for the association analysis of CSF Aβ42/Aβ40, Aβ PET and CSF p-Tau/Aβ40 in late amyloidosis stage  
## main part

graphics.off()
rm(list=ls())
library(stringr)
library(dplyr)
library(rprojroot)
library(ggseg)
library(ggpubr)
library(sjPlot) # For tab_model
library(doBy)
library(tidyverse)
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
nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta)
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Diag_AV1451_v1)

####################################################################################################
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         (Stage_CSF_PET_Abeta == "CSF+/PET+")|(Stage_CSF_PET_Abeta == "CSF-/PET-"))
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$Stage_CSF_PET_Abeta=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$Stage_CSF_PET_Abeta,
         levels = c("CSF-/PET-","CSF+/PET+"))
nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$Stage_CSF_PET_Abeta)
Col_Stage23_Ref = c("#4DAF4A","#E41A1C")



table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$Stage_CSF_PET_Abeta)

# ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ CSF_ABETA42.ABETA40_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
GLM_ENTORHINAL_FTP_Vs_CSF_Abeta_Stage23_Ref = 
  glm(ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        CSF_ABETA42.ABETA40_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_ENTORHINAL_FTP_Vs_CSF_Abeta_Stage23_Ref)
tab_model(GLM_ENTORHINAL_FTP_Vs_CSF_Abeta_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r_value = r$estimate
Sta_GLM = summary(GLM_ENTORHINAL_FTP_Vs_CSF_Abeta_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_ENTORHINAL_FTP_Vs_CSF_Abeta_Stage23_Ref = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = CSF_ABETA42.ABETA40_Closest_AV1451_v1, 
             y = ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(0.015,0.13,0.025),limits = c(0.015,0.13)) +
  scale_y_continuous(breaks = seq(0.7,2.3,0.2),limits = c(0.7,2.3)) +
  annotate("text", x = 0.0725, y = Inf,
           label = expression(bold(italic(beta)*' = -0.60, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "bold", color="black",vjust = 2) +
  xlab(expression('CSF A'*beta[42]*'/A'*beta[40])) +
  ylab("Entorhinal FTP SUVR") +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ PET_ABETA_CL_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$PET_ABETA_CL_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_ENTORHINAL_FTP_Vs_PET_Abeta_Stage23_Ref = 
  glm(ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        PET_ABETA_CL_Closest_AV1451_v1 +
        Gender + Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_ENTORHINAL_FTP_Vs_PET_Abeta_Stage23_Ref)
tab_model(GLM_ENTORHINAL_FTP_Vs_PET_Abeta_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$PET_ABETA_CL_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
r$conf.int
Sta_GLM = summary(GLM_ENTORHINAL_FTP_Vs_PET_Abeta_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_ENTORHINAL_FTP_Vs_PET_Abeta_Stage23_Ref = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = PET_ABETA_CL_Closest_AV1451_v1, 
             y = ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(-50,250,50),limits = c(-50,250)) +
  scale_y_continuous(breaks = seq(0.7,2.3,0.2),limits = c(0.7,2.3)) +
  annotate("text", x =  100, y = Inf,
           label = expression(bold(italic(beta)*' = 0.66, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "plain", color="black",vjust = 2) +
  xlab(expression('A'*beta*' PET (Centiloid)')) +
  ylab(" ") +
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank())



# ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ CSF_PTAU.ABETA40_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
GLM_ENTORHINAL_FTP_Vs_CSF_PTAU_Stage23_Ref = 
  glm(ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        CSF_PTAU.ABETA40_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_ENTORHINAL_FTP_Vs_CSF_PTAU_Stage23_Ref)
tab_model(GLM_ENTORHINAL_FTP_Vs_CSF_PTAU_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r_value = r$estimate
Sta_GLM = summary(GLM_ENTORHINAL_FTP_Vs_CSF_PTAU_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_ENTORHINAL_FTP_Vs_CSF_PTAU_Stage23_Ref = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = CSF_PTAU.ABETA40_Closest_AV1451_v1, 
             y = ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(0.0005,0.005,0.001),limits = c(0.0005,0.005)) +
  scale_y_continuous(breaks = seq(0.8,2.3,0.3),limits = c(0.8,2.3)) +
  annotate("text", x = 0.00275, y = Inf,
           label = expression(bold(italic(beta)*' = 0.73, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "plain", color="black",vjust = 2) +
  xlab(expression('CSF p-Tau'*'/A'*beta[40])) +
  ylab(" ") +
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank())

Plot_ENTORHINAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref = 
  ggarrange(
    Plot_ENTORHINAL_FTP_Vs_CSF_Abeta_Stage23_Ref,
    Plot_ENTORHINAL_FTP_Vs_PET_Abeta_Stage23_Ref,
    Plot_ENTORHINAL_FTP_Vs_CSF_PTAU_Stage23_Ref,
    nrow=1,ncol = 3, labels = c("A", "B","C"),widths = c(0.45,0.4),
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend="bottom")
# ggsave("Plot_ENTORHINAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref.jpg",
#        Plot_ENTORHINAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref, 
#        units="in", width=11, height=4, dpi=600)
# ggsave("Plot_ENTORHINAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref.pdf",
#        Plot_ENTORHINAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref, 
#        units="in", width=6.8, height=4, dpi=600)

# BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ CSF_ABETA42.ABETA40_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_INFERIORTEMPORAL_FTP_Vs_CSF_Abeta_Stage23_Ref = 
  glm(BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        CSF_ABETA42.ABETA40_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_INFERIORTEMPORAL_FTP_Vs_CSF_Abeta_Stage23_Ref)
tab_model(GLM_INFERIORTEMPORAL_FTP_Vs_CSF_Abeta_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
r$conf.int
Sta_GLM = summary(GLM_INFERIORTEMPORAL_FTP_Vs_CSF_Abeta_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_INFERIORTEMPORAL_FTP_Vs_CSF_Abeta_Stage23_Ref = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = CSF_ABETA42.ABETA40_Closest_AV1451_v1, 
             y = BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(0.015,0.13,0.025),limits = c(0.015,0.13)) +
  scale_y_continuous(breaks = seq(0.8,2.3,0.3),limits = c(0.8,2.3)) +
  annotate("text", x = 0.0725, y = Inf,
           label = expression(bold(italic(beta)*' = -0.48, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "plain", color="black",vjust = 2 ) +
  xlab(expression('CSF A'*beta[42]*'/A'*beta[40])) +
  labs(y=expression(paste(Braak[III/IV],' FTP SUVR')))+
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ PET_ABETA_CL_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$PET_ABETA_CL_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_INFERIORTEMPORAL_FTP_Vs_PET_Abeta_Stage23_Ref = 
  glm(BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        PET_ABETA_CL_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_INFERIORTEMPORAL_FTP_Vs_PET_Abeta_Stage23_Ref)
tab_model(GLM_INFERIORTEMPORAL_FTP_Vs_PET_Abeta_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$PET_ABETA_CL_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
r$conf.int
Sta_GLM = summary(GLM_INFERIORTEMPORAL_FTP_Vs_PET_Abeta_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_INFERIORTEMPORAL_FTP_Vs_PET_Abeta_Stage23_Ref_Non_Y_Label = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = PET_ABETA_CL_Closest_AV1451_v1, 
             y = BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(-50,250,50),limits = c(-50,250)) +
  scale_y_continuous(breaks = seq(0.8,2.3,0.3),limits = c(0.8,2.3)) +
  annotate("text", x =  100, y = Inf,
           label = expression(bold(italic(beta)*' = 0.54, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "bold", color="black",vjust = 2) +
  xlab(expression('A'*beta*' PET (Centiloid)')) +
  ylab(" ") +
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank())



# BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ CSF_PTAU.ABETA40_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_INFERIORTEMPORAL_FTP_Vs_PTAU_Stage23_Ref = 
  glm(BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        CSF_PTAU.ABETA40_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_INFERIORTEMPORAL_FTP_Vs_PTAU_Stage23_Ref)
tab_model(GLM_INFERIORTEMPORAL_FTP_Vs_PTAU_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
r$conf.int
Sta_GLM = summary(GLM_INFERIORTEMPORAL_FTP_Vs_PTAU_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_INFERIORTEMPORAL_FTP_Vs_PTAU_Stage23_Ref_Non_Y_Label = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = CSF_PTAU.ABETA40_Closest_AV1451_v1, 
             y = BRAAK34_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(0.0005,0.005,0.001),limits = c(0.0005,0.005)) +
  scale_y_continuous(breaks = seq(0.8,2.3,0.3),limits = c(0.8,2.3)) +
  annotate("text", x =  0.00275, y = Inf,
           label = expression(bold(italic(beta)*' = 0.75, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "plain", color="black",vjust = 2) +
  xlab(expression('CSF p-Tau'*'/A'*beta[40])) +
  ylab(" ") +
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank())


Plot_INFERIORTEMPORAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref = 
  ggarrange(
    Plot_INFERIORTEMPORAL_FTP_Vs_CSF_Abeta_Stage23_Ref,
    Plot_INFERIORTEMPORAL_FTP_Vs_PET_Abeta_Stage23_Ref_Non_Y_Label,
    Plot_INFERIORTEMPORAL_FTP_Vs_PTAU_Stage23_Ref_Non_Y_Label,
    nrow=1,ncol = 3,labels = c("D", "E","F"),widths = c(0.45,0.4),
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend="bottom")
# ggsave("Plot_INFERIORTEMPORAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref.jpg",
#        Plot_INFERIORTEMPORAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref, 
#        units="in", width=11, height=4, dpi=600)
# ggsave("Plot_INFERIORTEMPORAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref.pdf",
#        Plot_INFERIORTEMPORAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref, 
#        units="in", width=11, height=4, dpi=600)




# BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ CSF_ABETA42.ABETA40_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_SUPERIORFRONTAL_FTP_Vs_CSF_Abeta_Stage23_Ref = 
  glm(BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        CSF_ABETA42.ABETA40_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_SUPERIORFRONTAL_FTP_Vs_CSF_Abeta_Stage23_Ref)
tab_model(GLM_SUPERIORFRONTAL_FTP_Vs_CSF_Abeta_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_ABETA42.ABETA40_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
Sta_GLM = summary(GLM_SUPERIORFRONTAL_FTP_Vs_CSF_Abeta_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
Plot_SUPERIORFRONTAL_FTP_Vs_CSF_Abeta_Stage23_Ref = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = CSF_ABETA42.ABETA40_Closest_AV1451_v1, 
             y = BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(0.015,0.13,0.025),limits = c(0.015,0.13)) +
  scale_y_continuous(breaks = seq(0.8,2,0.2),limits = c(0.8,2)) +
  annotate("text", x = 0.0725, y = Inf,
           label = expression(bold(italic(beta)*' = -0.45, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "plain", color="black",vjust = 2) +
  xlab(expression('CSF A'*beta[42]*'/A'*beta[40])) +
  labs(y=expression(paste(Braak[V/VI],' FTP SUVR')))+
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ PET_ABETA_CL_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$PET_ABETA_CL_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_SUPERIORFRONTAL_FTP_Vs_PET_Abeta_Stage23_Ref = 
  glm(BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        PET_ABETA_CL_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_SUPERIORFRONTAL_FTP_Vs_PET_Abeta_Stage23_Ref)
tab_model(GLM_SUPERIORFRONTAL_FTP_Vs_PET_Abeta_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$PET_ABETA_CL_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
r$conf.int
Sta_GLM = summary(GLM_SUPERIORFRONTAL_FTP_Vs_PET_Abeta_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_SUPERIORFRONTAL_FTP_Vs_PET_Abeta_Stage23_Ref_Non_Y_Label = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = PET_ABETA_CL_Closest_AV1451_v1, 
             y = BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(-50,250,50),limits = c(-50,250)) +
  scale_y_continuous(breaks = seq(0.8,2,0.2),limits = c(0.8,2)) +
  annotate("text", x =  100, y = Inf,
           label = expression(bold(italic(beta)*' = 0.54, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "bold", color="black",vjust = 2) +
  xlab(expression('A'*beta*' PET (Centiloid)')) +
  ylab(" ") +
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank())



# BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 ~ CSF_PTAU.ABETA40_Closest_AV1451_v1
####################################################################################################
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1)
each(min,max)(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)
####################################################################################################
GLM_SUPERIORFRONTAL_FTP_Vs_PTAU_Stage23_Ref = 
  glm(BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1~
        CSF_PTAU.ABETA40_Closest_AV1451_v1+
        Gender+Age_AV1451_v1,
      family = gaussian(link = "identity"),
      data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref)
vif(GLM_SUPERIORFRONTAL_FTP_Vs_PTAU_Stage23_Ref)
tab_model(GLM_SUPERIORFRONTAL_FTP_Vs_PTAU_Stage23_Ref,
          show.se=TRUE,
          show.ci = 0.95,
          show.std = TRUE,
          string.std = "std. Beta",
          # string.se = "SE",
          digits=4,
          emph.p=TRUE,
          digits.p=5,
          string.p = "p value")
r = cor.test(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1, 
             Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref$CSF_PTAU.ABETA40_Closest_AV1451_v1, 
             method = "pearson", conf.level = 0.95)
r_value = r$estimate
r$conf.int
Sta_GLM = summary(GLM_SUPERIORFRONTAL_FTP_Vs_PTAU_Stage23_Ref)
p_value = Sta_GLM$coefficients[2,4]
r_p_value = paste(paste("R = ",as.character(signif(r_value,2)),sep = ""),
                  paste("p = ",as.character(signif(p_value,2)),sep = ""),sep = ", ")
r_p_value
Plot_SUPERIORFRONTAL_FTP_Vs_PTAU_Stage23_Ref_Non_Y_Label = 
  ggplot(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_Stage23_Ref,
         aes(x = CSF_PTAU.ABETA40_Closest_AV1451_v1, 
             y = BRAAK56_SUVR_FTP_Inf_CER_GM_Non_PVC_v1)) +
  geom_smooth(method="lm",  size = 1,linetype ="solid",
              color= "black", alpha = 0.5)+
  geom_point(aes(color= Stage_CSF_PET_Abeta,shape = Diag_AV1451_v1),alpha = 0.7) +
  scale_colour_manual(values = Col_Stage23_Ref)+
  scale_shape_manual(values = c(19,1)) +
  scale_x_continuous(breaks = seq(0.0005,0.005,0.001),limits = c(0.0005,0.005)) +
  scale_y_continuous(breaks = seq(0.8,2,0.2),limits = c(0.8,2)) +
  annotate("text", x =  0.00275, y = Inf,
           label = expression(bold(italic(beta)*' = 0.69, '*italic(p)*' < 0.001')),
           size = 3.5,fontface = "plain", color="black",vjust = 2) +
  xlab(expression('CSF p-Tau'*'/A'*beta[40])) +
  ylab(" ") +
  # guides(colour = guide_legend(nrow =2)) +
  theme(legend.title=element_blank())+
  theme(title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=10,face="plain"),
        axis.text=element_text(size=10),axis.title=element_text(size=10,face="plain"),
        legend.text=element_text(size=10),legend.position = "bottom", 
        legend.direction = "horizontal")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank())

Plot_SUPERIORFRONTAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref = 
  ggarrange(
    Plot_SUPERIORFRONTAL_FTP_Vs_CSF_Abeta_Stage23_Ref,
    Plot_SUPERIORFRONTAL_FTP_Vs_PET_Abeta_Stage23_Ref_Non_Y_Label,
    Plot_SUPERIORFRONTAL_FTP_Vs_PTAU_Stage23_Ref_Non_Y_Label,
    nrow=1,ncol = 3, labels = c("G", "H","I"),widths = c(0.45,0.4),
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend="bottom")
# ggsave("Plot_SUPERIORFRONTAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref.jpg",
#        Plot_SUPERIORFRONTAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref, 
#        units="in", width=11, height=4, dpi=600)


####################################################################################################
## 拼图+保存
####################################################################################################
Plot_Association_BL_FTP_SUVR_Stages_CSF_PET_Abeta = 
  ggarrange(
    Plot_ENTORHINAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref,
    Plot_INFERIORTEMPORAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref,
    Plot_SUPERIORFRONTAL_FTP_Vs_CSF_PET_Abeta_Stage23_Ref,
    nrow=3,ncol = 1, 
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend="bottom")

ggsave("D:/Project/guo_20210206/CSF_PET_tau/Output/main_results/with_otuliers/1_with_outliers/Plot_Association_4_new.tiff",
       Plot_Association_BL_FTP_SUVR_Stages_CSF_PET_Abeta,
       units="in", width=11, height=11, dpi=600)


