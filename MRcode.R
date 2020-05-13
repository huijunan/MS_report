# BL to posttx analysis

library(dplyr)
library(janitor)
library(ggplot2)
library(psycho)
#### Read data 
raw <- read.csv("data_bl_to_posttx.csv") # 270 patients at first
outcome <- read.csv("rise_widescales_dec18.csv") %>%
  select(ID, SF36V1_VITALITY_FUPTX, SF36V1_VITALITY_FU18mo) %>%
  rename(id = ID) # select only outcome and identifier

data <- raw %>%
  left_join(outcome, by = "id") %>% # put outcome variable in the dataset
  filter(treatment4 == 1 | 
           treatment4 == 2 | 
           treatment4 == 3) %>% # selecting only chemo/radiation patients (206)
  filter(!is.na(SF36V1_VITALITY_FUPTX)) # opt out patient without outcome (194)

immune <- read.csv("rise_immune_bl.csv") %>%
  filter(id %in% data$id) %>% # select patients only in dataset
  select(-time, -IMMUNE_DT,-deleted_ifn_g_value, 
         -deleted_il_10_value, -deleted_il_6_value, -deleted_il_8_value, 
         -deleted_tnf_a_value, -deleted_crp_value, -deleted_stnfr2_value)
immune_std_temp <- immune %>% 
  select(id, LOG_IMMUNE_IFN_G, LOG_IMMUNE_IL_10, LOG_IMMUNE_IL_6, LOG_IMMUNE_IL_8,
           LOG_IMMUNE_TNF_A, LOG_IMMUNE_CRP, LOG_IMMUNE_sTNFR2) %>%
  mutate_at(scale, .vars = vars(-id))
immune_std_temp$immune_tot <- rowSums(immune_std_temp %>% select(-id))
immune_std <- immune_std_temp %>% select(id, immune_tot)


snp <- read.csv("immune_long_analytic.csv") %>%
  filter(time == "0 BASELINE") %>%
  select(id, TNF_highalleles, IL6_highalleles, IL1B_highalleles) %>%
  distinct()


### Check patient uniqueness
# data %>% group_by(id) %>% summarize(n_distinct(id)) %>% View()
# immune %>% group_by(id) %>% summarize(n_distinct(id)) %>% View()
# outcome %>% group_by(ID) %>% summarize(n_distinct(ID)) %>% View()
# 
# immune %>% filter(!id %in% data$id) %>% View()
# outcome %>% filter(!ID %in% data$id) %>% View()
# data %>% filter(!id %in% outcome$ID) %>% View()

#### Data cleaning
data1 <- data %>%
  left_join(snp, by = "id") %>%
  mutate(Outcome = ifelse(SF36V1_VITALITY_FUPTX <= 45, 1, 0),
         Vital_BL = ifelse(SF36V1_VITALITY_BL <= 45, 1, 0)) %>%
  # mutate(stage = case_when(.$STAGEDX_PG == 2 ~ 1, 
  #                          .$STAGEDX_PG == 3 ~ 1, 
  #                          .$STAGEDX_PG == 1 ~ 0,
  #                          .$STAGEDX_PG == 0 ~ 0)) %>%
  select(id, SF36V1_VITALITY_BL, Vital_BL, MFSI_GEN_BL, AGE_BL, race5, educ3cat,
         IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, SPS_ATTACH_BL,
         employedYN, MARRIED2, SURGTYPE_ENROLLMENT, treatment4, BMI_BL, charlson, ctq3cat, 
         SCID_PHMDD, TNF_highalleles, IL6_highalleles, IL1B_highalleles, Outcome) %>%
  # make dummy variables
  mutate(white = ifelse(race5 == 0, 1, 0),
         edulow = ifelse(educ3cat == 1, 1, 0),
         edumed = ifelse(educ3cat == 2, 1, 0),
         chemo = ifelse(treatment4 == 1, 1, 0),
         radia = ifelse(treatment4 == 2, 1, 0),
         # char1 = ifelse(charlson == 1, 1, 0),
         # char23 = ifelse(charlson == 2 | charlson == 3, 1, 0),
         ctq = ifelse(ctq3cat == 0, 0, 1),
         SPS_max = ifelse(SPS_ATTACH_BL == 16, 1, 0),
         tnf2cat = ifelse(TNF_highalleles == 2, 1, 0),
         il6_2cat = ifelse(IL6_highalleles == 2, 1, 0),
         il1b_2cat = ifelse(IL6_highalleles == 0, 0, 1),
         lump = ifelse(SURGTYPE_ENROLLMENT == 1, 1, 0),
         mast = ifelse(SURGTYPE_ENROLLMENT == 2 | SURGTYPE_ENROLLMENT == 3 | SURGTYPE_ENROLLMENT == 4, 1, 0)) %>% 
  left_join(immune, by = "id") %>%
  left_join(immune_std, by = "id")

immune1 <- immune %>%
  left_join(data1 %>% select(id, Outcome), by = "id")

fat <- data1 %>% filter(Outcome == 1)
nonfat <- data1 %>% filter(Outcome == 0)

fat_immune <- immune1 %>% filter(Outcome == 1)
nonfat_immune <- immune1 %>% filter(Outcome == 0)

#### Descriptive Analysis
# continuous: including age
mapply(mean, data1 %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                              IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson), na.rm=T)
mapply(mean, fat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, IES_BL, 
                            FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson), na.rm=T)
mapply(mean, nonfat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                               IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson), na.rm=T)

mapply(sd, data1 %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                            IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson), na.rm=T)
mapply(sd, fat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                          IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson), na.rm=T)
mapply(sd, nonfat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                             IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson), na.rm=T)

mapply(mean, immune %>% select(-id))
mapply(mean, fat_immune %>% select(-id))
mapply(mean, nonfat_immune %>% select(-id))

mapply(sd, immune %>% select(-id))
mapply(sd, fat_immune %>% select(-id))
mapply(sd, nonfat_immune %>% select(-id))


# categorical: including Vital_BL, 
# race5, income3, educ3cat,
# employedYN, MARRIED2, stage, surg_post, treatment4, 
# charlson, ctq3cat, 
# SCID_PHMDD, SNP_IL6_LETTER, SNP_TNF_LETTER, SNP_IL1B_LETTER
mod1 <- glm(Outcome ~ 1, family = "binomial", data = data1) # with rank
tabyl(data1$Vital_BL)
tabyl(fat$Vital_BL)
tabyl(nonfat$Vital_BL)

Vital_BL <- glm(formula = Outcome ~ Vital_BL, family = "binomial", data = data1)
summary(Vital_BL)
#------------
tabyl(data1$white)
tabyl(fat$white)
tabyl(nonfat$white)

race <- glm(formula = Outcome ~ white, family = "binomial", data = data1)
summary(race)
#------------
# tabyl(data1$income3)
# tabyl(fat$income3)
# tabyl(nonfat$income3)
# income <- glm(formula = Outcome ~ incomelow + incomemed, family = "binomial", data = data1)
# summary(income)
# mod2 <- glm(Outcome ~ 1, family = "binomial", data = data1 %>% filter(!is.na(income3))) # with rank
# anova(income, mod2, test="LRT")
#------------
bmi_cat <- glm(formula = Outcome ~ bmi_low + bmi_ow + bmi_ob, family = "binomial", data = data1)
summary(bmi_cat)
anova(bmi_cat, mod1, test="LRT")
#------------
tabyl(data1$educ3cat)
tabyl(fat$educ3cat)
tabyl(nonfat$educ3cat)
edu <- glm(formula = Outcome ~ edulow + edumed, family = "binomial", data = data1)
summary(edu)
anova(edu, mod1, test="LRT")
#------------
tabyl(data1$employedYN)
tabyl(fat$employedYN)
tabyl(nonfat$employedYN)
employedYN <- glm(formula = Outcome ~ employedYN, family = "binomial", data = data1)
summary(employedYN)
#------------
tabyl(data1$MARRIED2)
tabyl(fat$MARRIED2)
tabyl(nonfat$MARRIED2)
MARRIED2 <- glm(formula = Outcome ~ MARRIED2, family = "binomial", data = data1)
summary(MARRIED2)
#------------
tabyl(data1$tnf2cat)
tabyl(fat$tnf2cat)
tabyl(nonfat$tnf2cat)
tnf2cat <- glm(formula = Outcome ~ tnf2cat, family = "binomial", data = data1)
summary(tnf2cat)
#------------
tabyl(data1$il6_2cat)
tabyl(fat$il6_2cat)
tabyl(nonfat$il6_2cat)
il6_2cat <- glm(formula = Outcome ~ il6_2cat, family = "binomial", data = data1)
summary(il6_2cat)
#------------
tabyl(data1$il1b_2cat)
tabyl(fat$il1b_2cat)
tabyl(nonfat$il1b_2cat)
il1b_2cat <- glm(formula = Outcome ~ il1b_2cat, family = "binomial", data = data1)
summary(il1b_2cat)
#------------
tabyl(data1$SPS_max)
tabyl(fat$SPS_max)
tabyl(nonfat$SPS_max)
SPS_max <- glm(formula = Outcome ~ SPS_max, family = "binomial", data = data1)
summary(SPS_max)
#------------
# tabyl(data1$stage)
# tabyl(fat$stage)
# tabyl(nonfat$stage)
# stage <- glm(formula = Outcome ~ stage, family = "binomial", data = data1)
# summary(stage)
#------------
tabyl(data1$SURGTYPE_ENROLLMENT)
tabyl(fat$SURGTYPE_ENROLLMENT)
tabyl(nonfat$surg_post)
surg <- glm(formula = Outcome ~ lump + mast, family = "binomial", data = data1)
summary(surg)
anova(surg, mod1, test="LRT")
#------------
tabyl(data1$treatment4)
tabyl(fat$treatment4)
tabyl(nonfat$treatment4)
treatment <- glm(formula = Outcome ~ chemo + radia, family = "binomial", data = data1)
summary(treatment)
anova(treatment, mod1, test="LRT")
#------------
tabyl(data1$charlson)
tabyl(fat$charlson)
tabyl(nonfat$charlson)
charlson <- glm(formula = Outcome ~ charlson, family = "binomial", data = data1)
summary(charlson)
#------------
tabyl(data1$ctq)
tabyl(fat$ctq)
tabyl(nonfat$ctq)
ctq <- glm(formula = Outcome ~ ctq, family = "binomial", data = data1)
summary(ctq)
#------------
tabyl(data1$SCID_PHMDD)
tabyl(fat$SCID_PHMDD)
tabyl(nonfat$SCID_PHMDD)
SCID_PHMDD <- glm(formula = Outcome ~ SCID_PHMDD, family = "binomial", data = data1)
summary(SCID_PHMDD)
#------------
tabyl(data1$SNP_TNF_LETTER)
tabyl(fat$SNP_TNF_LETTER)
tabyl(nonfat$SNP_TNF_LETTER)
TNF <- glm(formula = Outcome ~ TNF_AA + TNF_GA, family = "binomial", data = data1)
summary(TNF)
anova(TNF, mod1, test="LRT")
#------------
tabyl(data1$SNP_IL6_LETTER)
tabyl(fat$SNP_IL6_LETTER)
tabyl(nonfat$SNP_IL6_LETTER)
IL6 <- glm(formula = Outcome ~ IL6_CC + IL6_GC, family = "binomial", data = data1)
summary(IL6)
anova(IL6, mod1, test="LRT")
#------------
tabyl(data1$SNP_IL1B_LETTER)
tabyl(fat$SNP_IL1B_LETTER)
tabyl(nonfat$SNP_IL1B_LETTER)
IL1B <- glm(Outcome ~ IL1B_AA + IL1B_AG, family = "binomial", data = data1)
summary(IL1B)
anova(IL1B, mod1, test="LRT")
#------------


# simple logistic regressions
IES <- glm(formula = Outcome ~ IES_BL, family = "binomial", data = data1)
summary(IES)

FCS <- glm(formula = Outcome ~ FCS_BL, family = "binomial", data = data1)
summary(FCS)

PSS <- glm(formula = Outcome ~ PSS_BL, family = "binomial", data = data1)
summary(PSS)

PSQI <- glm(formula = Outcome ~ PSQI_BL, family = "binomial", data = data1)
summary(PSQI)

CESD <- glm(formula = Outcome ~ CESD_BL, family = "binomial", data = data1)
summary(CESD)

SPS <- glm(formula = Outcome ~ SPS_ATTACH_BL, family = "binomial", data = data1)
summary(SPS)

age <- glm(formula = Outcome ~ AGE_BL, family = "binomial", data = data1)
summary(age)

sfv_bl <- glm(formula = Outcome ~ SF36V1_VITALITY_BL, family = "binomial", data = data1)
summary(sfv_bl)

bmi <- glm(formula = Outcome ~ BMI_BL, family = "binomial", data = data1)
summary(bmi)

mfsi_bl <- glm(formula = Outcome ~ MFSI_GEN_BL, family = "binomial", data = data1)
summary(mfsi_bl)

LOG_IMMUNE_IFN_G <- glm(formula = Outcome ~ LOG_IMMUNE_IFN_G, family = "binomial", data = data1)
summary(LOG_IMMUNE_IFN_G)

LOG_IMMUNE_IL_10 <- glm(formula = Outcome ~ LOG_IMMUNE_IL_10, family = "binomial", data = data1)
summary(LOG_IMMUNE_IL_10)

LOG_IMMUNE_IL_6 <- glm(formula = Outcome ~ LOG_IMMUNE_IL_6, family = "binomial", data = data1)
summary(LOG_IMMUNE_IL_6)

LOG_IMMUNE_IL_8 <- glm(formula = Outcome ~ LOG_IMMUNE_IL_8, family = "binomial", data = data1)
summary(LOG_IMMUNE_IL_8)

LOG_IMMUNE_TNF_A <- glm(formula = Outcome ~ LOG_IMMUNE_TNF_A, family = "binomial", data = data1)
summary(LOG_IMMUNE_TNF_A)

LOG_IMMUNE_CRP <- glm(formula = Outcome ~ LOG_IMMUNE_CRP, family = "binomial", data = data1)
summary(LOG_IMMUNE_CRP)

LOG_IMMUNE_sTNFR2 <- glm(formula = Outcome ~ LOG_IMMUNE_sTNFR2, family = "binomial", data = data1)
summary(LOG_IMMUNE_sTNFR2)

IMMUNE_IFN_G <- glm(formula = Outcome ~ IMMUNE_IFN_G, family = "binomial", data = data1)
summary(IMMUNE_IFN_G)

IMMUNE_IL_10 <- glm(formula = Outcome ~ IMMUNE_IL_10, family = "binomial", data = data1)
summary(IMMUNE_IL_10)

IMMUNE_IL_6 <- glm(formula = Outcome ~ IMMUNE_IL_6, family = "binomial", data = data1)
summary(IMMUNE_IL_6)

IMMUNE_IL_8 <- glm(formula = Outcome ~ IMMUNE_IL_8, family = "binomial", data = data1)
summary(IMMUNE_IL_8)

IMMUNE_TNF_A <- glm(formula = Outcome ~ IMMUNE_TNF_A, family = "binomial", data = data1)
summary(IMMUNE_TNF_A)

IMMUNE_CRP <- glm(formula = Outcome ~ IMMUNE_CRP, family = "binomial", data = data1)
summary(IMMUNE_CRP)

IMMUNE_sTNFR2 <- glm(formula = Outcome ~ IMMUNE_sTNFR2, family = "binomial", data = data1)
summary(IMMUNE_sTNFR2)

immune_tot <- glm(formula = Outcome ~ immune_tot, family = "binomial", data = data1)
summary(immune_tot)


# log-transformed immune markers distributions
ggplot(data1, aes(x=LOG_IMMUNE_IFN_G))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_IFN_G") 

ggplot(data1, aes(x=LOG_IMMUNE_IL_10))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_IL_10") 

ggplot(data1, aes(x=LOG_IMMUNE_IL_6))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_IL_6") 

ggplot(data1, aes(x=LOG_IMMUNE_IL_8))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_IL_8") 

ggplot(data1, aes(x=LOG_IMMUNE_TNF_A))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_TNF_A")

ggplot(data1, aes(x=LOG_IMMUNE_CRP))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_CRP")

ggplot(data1, aes(x=LOG_IMMUNE_sTNFR2))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-5.5:9)) + 
  ggtitle("Distribution of LOG_IMMUNE_sTNFR2")

data.glm <- na.omit(data1)
# overall logistic regression
glm <- glm(formula = Outcome ~ AGE_BL + BMI_BL + Vital_BL + white  + 
             edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
             ctq + SCID_PHMDD + 
             IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
             tnf2cat + il6_2cat + il1b_2cat +
             LOG_IMMUNE_IFN_G + LOG_IMMUNE_IL_10 + LOG_IMMUNE_IL_6 + LOG_IMMUNE_IL_8 +
             LOG_IMMUNE_TNF_A + LOG_IMMUNE_CRP + LOG_IMMUNE_sTNFR2, family = "binomial", data = data.glm)
summary(glm)
glm2 <- glm(formula = Outcome ~ AGE_BL + BMI_BL + Vital_BL + white  + 
             edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
              ctq + SCID_PHMDD + 
              IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
             tnf2cat + il6_2cat + il1b_2cat +
             immune_tot, family = "binomial", data = data.glm)
summary(glm2) # using immune score
glm3 <- glm(formula = Outcome ~ AGE_BL + BMI_BL + I(BMI_BL^2) + Vital_BL + white  + 
              edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
              ctq + SCID_PHMDD + 
              IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
              tnf2cat + il6_2cat + il1b_2cat +
              immune_tot, family = "binomial", data = data.glm)
summary(glm3) # using quadratic bmi

# mlr with each immune marker at a time
IFN_G <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
               edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
               IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
               ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
               BMI_BL + LOG_IMMUNE_IFN_G, family = "binomial", data = data1)
summary(IFN_G)
IL_10 <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
               edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
               IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
               ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
               BMI_BL + LOG_IMMUNE_IL_10, family = "binomial", data = data1)
summary(IL_10)

IL_6 <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
              edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
              IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
              ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
              BMI_BL + LOG_IMMUNE_IL_6, family = "binomial", data = data1)
summary(IL_6)
IL_8 <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
              edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
              IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
              ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
              BMI_BL + LOG_IMMUNE_IL_8, family = "binomial", data = data1)
summary(IL_8)
TNF_A <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
               edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
               IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
               ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
               BMI_BL + LOG_IMMUNE_TNF_A, family = "binomial", data = data1)
summary(TNF_A)
IMMUNE_CRP <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
                    edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
                    IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
                    ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
                    BMI_BL + LOG_IMMUNE_CRP, family = "binomial", data = data1)
summary(IMMUNE_CRP)
sTNFR2 <- glm(formula = Outcome ~ AGE_BL + Vital_BL + white  + 
                edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
                IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
                ctq + SCID_PHMDD + tnf2cat + il6_2cat + il1b_2cat +
                BMI_BL + LOG_IMMUNE_sTNFR2, family = "binomial", data = data1)
summary(sTNFR2)


# immune markers correlation table
round(cor(immune %>% select(LOG_IMMUNE_IFN_G, LOG_IMMUNE_IL_10, LOG_IMMUNE_IL_6, LOG_IMMUNE_IL_8,
                        LOG_IMMUNE_TNF_A, LOG_IMMUNE_CRP, LOG_IMMUNE_sTNFR2)),2) %>% View()
# PCA
pca <- prcomp(immune %>% select(LOG_IMMUNE_IFN_G, LOG_IMMUNE_IL_10, LOG_IMMUNE_IL_6, LOG_IMMUNE_IL_8,
                                LOG_IMMUNE_TNF_A, LOG_IMMUNE_CRP, LOG_IMMUNE_sTNFR2), center = TRUE,scale. = TRUE)
summary(pca)
# plot pca
data2 <- data1 %>% filter(!is.na(LOG_IMMUNE_CRP)) %>% mutate(Outcome = as.factor(Outcome))
library("factoextra")
fviz_pca_ind(pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = data2$Outcome, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Fatigue") +
  ggtitle("2D PCA-plot from 7 log transformed immune markers") +
  theme(plot.title = element_text(hjust = 0.5))

# ANOVA for SNPs
data4 <- data1 %>% filter(il1b_2cat != "NA") %>% filter(!is.na(LOG_IMMUNE_IL_6))
res.aov1 <- aov(LOG_IMMUNE_TNF_A ~ tnf2cat, data = data4)
summary(res.aov1)
ggplot(data4, aes(x=tnf2cat, y=LOG_IMMUNE_TNF_A, fill=tnf2cat))+
       geom_boxplot() + ggtitle("Plot of Log_TNF_A by SNP_TNF, ANOVA p-value = 0.47")

res.aov2 <- aov(LOG_IMMUNE_IL_6 ~ il6_2cat, data = data4)
summary(res.aov2)

res.aov3 <- aov(LOG_IMMUNE_CRP ~ il6_2cat, data = data4)
summary(res.aov3)

summary(aov(LOG_IMMUNE_IFN_G ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_10 ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_6 ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_8 ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_TNF_A ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_CRP ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_sTNFR2 ~ tnf2cat, data = data4))
summary(aov(immune_tot ~ tnf2cat, data = data4))
        
summary(aov(LOG_IMMUNE_IFN_G ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_10 ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_6 ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_8 ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_TNF_A ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_CRP ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_sTNFR2 ~ il6_2cat, data = data4))
summary(aov(immune_tot ~ il6_2cat, data = data4))

summary(aov(LOG_IMMUNE_IFN_G ~ il1b_2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_10 ~ il1b_2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_6 ~ il1b_2cat, data = data4))
summary(aov(LOG_IMMUNE_IL_8 ~ il1b_2cat, data = data4))
summary(aov(LOG_IMMUNE_TNF_A ~ il1b_2cat, data = data4))
summary(aov(LOG_IMMUNE_CRP ~ il1b_2cat, data = data4))
summary(aov(LOG_IMMUNE_sTNFR2 ~ il1b_2cat, data = data4))
summary(aov(immune_tot ~ il1b_2cat, data = data4))

# Variable selection

step1 <- stepAIC(glm, scope=list(upper = ~ AGE_BL + BMI_BL + Vital_BL + white  + 
                                   edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
                                   ctq + SCID_PHMDD + 
                                   IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
                                   tnf2cat + il6_2cat + il1b_2cat +
                                   LOG_IMMUNE_IFN_G + LOG_IMMUNE_IL_10 + LOG_IMMUNE_IL_6 + LOG_IMMUNE_IL_8 +
                                  LOG_IMMUNE_TNF_A + LOG_IMMUNE_CRP + LOG_IMMUNE_sTNFR2, 
                                lower = ~ AGE_BL + BMI_BL + Vital_BL + white  + 
                                  edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
                                  ctq + SCID_PHMDD + 
                                  IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
                                  tnf2cat + il6_2cat + il1b_2cat), trace=FALSE)
summary(step1)

step2 <- stepAIC(glm)
summary(step2)
step2_f <- glm(formula = Outcome ~ BMI_BL + Vital_BL + 
                edulow + edumed + employedYN + charlson +
                PSQI_BL + CESD_BL + tnf2cat + LOG_IMMUNE_IL_6 + LOG_IMMUNE_CRP, family = "binomial", data = data1)
summary(step2_f)
