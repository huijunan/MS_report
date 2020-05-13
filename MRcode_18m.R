# BL to posttx(18months) analysis

library(dplyr)
library(janitor)
library(ggplot2)
library(psycho)

#### Read data 
raw <- read.csv("data_bl_to_posttx.csv") # 270 patients at first
outcome <- read.csv("rise_widescales_dec18.csv") %>%
  select(ID, SF36V1_VITALITY_FU18mo) %>%
  rename(id = ID) # select only outcome and identifier

data <- raw %>%
  left_join(outcome, by = "id") %>% # put outcome variable in the dataset
  filter(treatment4 == 1 | 
           treatment4 == 2 | 
           treatment4 == 3) %>% # selecting only chemo/radiation patients (206)
  filter(!is.na(SF36V1_VITALITY_FU18mo)) # opt out patient without outcome (186)

immune <- read.csv("rise_immune_bl.csv") %>%
  filter(id %in% data$id) %>% # select patients only in dataset
  select(id, LOG_IMMUNE_CRP)

snp <- read.csv("immune_long_analytic.csv") %>%
  filter(time == "0 BASELINE") %>%
  select(id, TNF_highalleles, IL6_highalleles, IL1B_highalleles) %>%
  distinct()

#### Data cleaning
data1 <- data %>%
  left_join(snp, by = "id") %>%
  mutate(Outcome = ifelse(SF36V1_VITALITY_FU18mo <= 45, 1, 0),
         Vital_BL = ifelse(SF36V1_VITALITY_BL <= 45, 1, 0)) %>%
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
         ctq = ifelse(ctq3cat == 0, 0, 1),
         SPS_max = ifelse(SPS_ATTACH_BL == 16, 1, 0),
         tnf2cat = ifelse(TNF_highalleles == 2, 1, 0),
         il6_2cat = ifelse(IL6_highalleles == 2, 1, 0),
         il1b_2cat = ifelse(IL6_highalleles == 0, 0, 1),
         lump = ifelse(SURGTYPE_ENROLLMENT == "1 Lumpectomy", 1, 0),
         mast = ifelse(SURGTYPE_ENROLLMENT == "2 Mastectomy WITHOUT immediate reconstruction" | 
                         SURGTYPE_ENROLLMENT == "3 Mastectomy WITH immediate reconstruction" | 
                         SURGTYPE_ENROLLMENT == "4 Delayed reconstruction of previous mastectomy", 1, 0)) %>% 
  left_join(immune, by = "id")

immune1 <- immune %>%
  left_join(data1 %>% select(id, Outcome), by = "id")

fat <- data1 %>% filter(Outcome == 1)
nonfat <- data1 %>% filter(Outcome == 0)

ggplot(data1, aes(x=LOG_IMMUNE_CRP))+
  geom_histogram(binwidth=0.25,color="black", fill="darkseagreen3",na.rm=T)+
  scale_x_continuous(breaks=c(-2.5:5)) + 
  ggtitle("Distribution of LOG_IMMUNE_CRP")

#### Descriptive Analysis
# continuous: including age

mapply(mean, data1 %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                              IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson, LOG_IMMUNE_CRP), na.rm=T)
mapply(mean, fat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, IES_BL, 
                            FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson, LOG_IMMUNE_CRP), na.rm=T)
mapply(mean, nonfat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                               IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson, LOG_IMMUNE_CRP), na.rm=T)

mapply(sd, data1 %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                            IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson, LOG_IMMUNE_CRP), na.rm=T)
mapply(sd, fat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                          IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson, LOG_IMMUNE_CRP), na.rm=T)
mapply(sd, nonfat %>% select(AGE_BL, BMI_BL, SF36V1_VITALITY_BL, MFSI_GEN_BL, 
                             IES_BL, FCS_BL, PSS_BL, PSQI_BL, CESD_BL, charlson, LOG_IMMUNE_CRP), na.rm=T)
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
tabyl(data1$SURGTYPE_ENROLLMENT)
tabyl(fat$SURGTYPE_ENROLLMENT)
tabyl(nonfat$SURGTYPE_ENROLLMENT)
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

age <- glm(formula = Outcome ~ AGE_BL, family = "binomial", data = data1)
summary(age)

sfv_bl <- glm(formula = Outcome ~ SF36V1_VITALITY_BL, family = "binomial", data = data1)
summary(sfv_bl)

bmi <- glm(formula = Outcome ~ BMI_BL, family = "binomial", data = data1)
summary(bmi)

mfsi_bl <- glm(formula = Outcome ~ MFSI_GEN_BL, family = "binomial", data = data1)
summary(mfsi_bl)

LOG_IMMUNE_CRP <- glm(formula = Outcome ~ LOG_IMMUNE_CRP, family = "binomial", data = data1)
summary(LOG_IMMUNE_CRP)

charlson <- glm(formula = Outcome ~ charlson, family = "binomial", data = data1)
summary(charlson)


data.glm <- na.omit(data1)
# overall logistic regression
glm <- glm(formula = Outcome ~ AGE_BL + BMI_BL + Vital_BL + white  + 
             edulow + edumed + employedYN + MARRIED2 + lump + mast + chemo + radia + charlson +
             ctq + SCID_PHMDD + 
             IES_BL + FCS_BL + PSS_BL + PSQI_BL + CESD_BL + SPS_max +
             tnf2cat + il6_2cat + il1b_2cat + LOG_IMMUNE_CRP, family = "binomial", data = data.glm)
summary(glm)

# ANOVA
data4 <- data1 %>% filter(il1b_2cat != "NA") %>% filter(!is.na(LOG_IMMUNE_CRP))
summary(aov(LOG_IMMUNE_CRP ~ tnf2cat, data = data4))
summary(aov(LOG_IMMUNE_CRP ~ il6_2cat, data = data4))
summary(aov(LOG_IMMUNE_CRP ~ il1b_2cat, data = data4))

# Variable selection
library(MASS)
step <- stepAIC(glm)
summary(step)
