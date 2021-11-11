graphics.off()
rm(list=ls())
library(readxl)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggseg)
library(ggpubr)
library(sjPlot)
library(plyr)
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
#####################################################################################################################
## Still cutoffs form LONI and GMM method for ABETA42/40 ratio
Cutoff_CSF_MASS_ABETA42 = 1079
Cutoff_CSF_MASS_ABETA42.ABETA40 = 0.138
Cutoff_CSF_MASS_ABETA42.ABETA38 = 0.627
#####################################################################################################################
Cutoff_SUVR_CER = 1.11
Cutoff_SUVR_Big_Ref = 0.82
Cutoff_SUVR_FBB_CER = 1.08
Cutoff_Centiloid = 24.4
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

# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A =
#   subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
#          !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 1.6)&(Stage_CSF_PET_Abeta == "CSF+/PET-")))
# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 =
#   subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A,
#          !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 2)&(Stage_CSF_PET_Abeta == "CSF+/PET+")))


Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1[is.na(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$
                                                         APOE_status),c("APOE_status")] = "APOE4_Non_carrier"
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta)

Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Gender=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Gender,
         levels = c("Female","Male"))
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$APOE_status=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$APOE_status,
         levels = c("APOE4_carrier","APOE4_Non_carrier"))
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Diag_AV1451_v1=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Diag_AV1451_v1,
         levels = c("CU","MCI"))
# table1(~ Diag_AV1451_v1+Age_AV1451_v1 + Gender + Education + 
#          APOE_status+aHCV_Closest_AV1451_v1+PACC_Digit_LONI_Closest_AV1451_v1 + 
#          CSF_ABETA42.ABETA40_Closest_AV1451_v1 + PET_ABETA_CL_Closest_AV1451_v1
#        |Stage_CSF_PET_Abeta,
#        data=Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)

library(compareGroups)
library(tcltk2)


# 
# Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 = rename(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
#             MCI = Diag_AV1451_v1,Age = Age_AV1451_v1,
#             Sex = Gender,APOE = APOE_status,aHCV = aHCV_Closest_AV1451_v1,
#             PACC = PACC_Digit_LONI_Closest_AV1451_v1,
#             CSF = CSF_ABETA42.ABETA40_Closest_AV1451_v1,
#             PET = PET_ABETA_CL_Closest_AV1451_v1)

restab <- descrTable(Stage_CSF_PET_Abeta~Diag_AV1451_v1+Age_AV1451_v1+Gender+APOE_status+Education+
                       aHCV_Closest_AV1451_v1+PACC_Digit_LONI_Closest_AV1451_v1 + 
                       CSF_ABETA42.ABETA40_Closest_AV1451_v1 + PET_ABETA_CL_Closest_AV1451_v1 +
                       CSF_PTAU.ABETA40_Closest_AV1451_v1,digits = 5,
            data = Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
            hide = c(Diag_AV1451_v1 = "CU",Gender = "Male",APOE_status = "APOE4_Non_carrier"),
            method = c(Age_AV1451_v1 = 2, Education = 2,aHCV_Closest_AV1451_v1 = 2,
                       PACC_Digit_LONI_Closest_AV1451_v1 = 2,
                       CSF_ABETA42.ABETA40_Closest_AV1451_v1 = 2,
                       PET_ABETA_CL_Closest_AV1451_v1 = 2,CSF_PTAU.ABETA40_Closest_AV1451_v1 = 2),show.p.mul = T)
export2html(restab, file='D:/Project/guo_20210206/CSF_PET_tau/Script/SUVR_defined_groups/table2.html')


## wilcox.test 的 estimate CI
## fisher exact test 的 odds ratio CI

descrTable(formula, data, subset, na.action = NULL, y = NULL, Xext = NULL, 
           selec = NA, method = 1, timemax = NA, alpha = 0.05, min.dis = 5, max.ylev = 5, 
           max.xlev = 10, include.label = TRUE, Q1 = 0.25, Q3 = 0.75, simplify = TRUE, 
           ref = 1, ref.no = NA, fact.ratio = 1, ref.y = 1, p.corrected = TRUE, 
           compute.ratio = TRUE, include.miss = FALSE, oddsratio.method = "midp", 
           chisq.test.perm = FALSE, byrow = FALSE, chisq.test.B = 2000, chisq.test.seed = NULL,
           Date.format = "d-mon-Y", var.equal = TRUE, conf.level = 0.95, surv = FALSE,
           riskratio = FALSE, riskratio.method = "wald", compute.prop = FALSE, 
           lab.missing = "'Missing'",
           hide = NA, digits = NA, type = NA, show.p.overall = TRUE,
           show.all, show.p.trend, show.p.mul = FALSE, show.n, show.ratio =
             FALSE, show.descr = TRUE, show.ci = FALSE, hide.no = NA, digits.ratio = NA,
           show.p.ratio = show.ratio, digits.p = 3, sd.type = 1, q.type = c(1, 1),
           extra.labels = NA, all.last = FALSE)

























