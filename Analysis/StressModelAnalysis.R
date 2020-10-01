# Stress Test Stat Models
# Set working directory and load necessary packages
library(effects)
library(lme4)
library(car)
library(MASS)
library(lmerTest)
library(multcomp)
library(MuMIn)

# Import data
StressMaster <- read.csv("Data/Stress_Master.csv")
StressMaster$group<-paste(StressMaster$Species, StressMaster$Origin, StressMaster$Transplant, StressMaster$Treatment)
set.seed(3028)

# PAM
# Eliminate missing PAM data
PAMData <- StressMaster[complete.cases(StressMaster$Tf), ] 
PAMData$Tf <- as.numeric(as.character(PAMData$Tf))

# Model selection
modelPAM1 <- lmer(Tf ~ Origin * Transplant * Species * Treatment + (1|Parent), data = PAMData, REML = F)
modelPAM2 <- lmer(log(Tf) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = PAMData, REML = F)
modelPAM3 <- lmer(sqrt(Tf) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = PAMData, REML = F)

model.sel(modelPAM1, modelPAM2, modelPAM3)
summary(modelPAM2)
anova(modelPAM2)
shapiro.test(residuals(modelPAM2))
qqPlot(residuals(modelPAM2))

# Gross Photosynthesis
# Eliminate missing GP data
GPData <- StressMaster[complete.cases(StressMaster$tfGPSA),]

# Model selection
GPSA1 = lmer(tfGPSA ~ Origin * Transplant * Species * Treatment + (1|Parent), data = GPData)
GPSA2 = lmer(log(tfGPSA) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = GPData)
GPSA3 = lmer(sqrt(tfGPSA) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = GPData)

model.sel(GPSA1, GPSA2, GPSA3)
summary(GPSA1)
anova(GPSA1)
shapiro.test(residuals(GPSA1))
qqPlot(residuals(GPSA1))

# Net Photosynthesis
# Eliminate missing NP data
NPData <- StressMaster[complete.cases(StressMaster$tfNPSA),]

# Model selection
NPSA1 = lmer(tfNPSA ~ Origin * Transplant * Species * Treatment + (1|Parent), data = NPData)
NPSA2 = lmer(log(tfNPSA) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = NPData)
NPSA3 = lmer(sqrt(tfNPSA) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = NPData)

model.sel(NPSA1, NPSA2, NPSA3)
summary(NPSA1)
anova(NPSA1)
shapiro.test(residuals(NPSA1))
qqPlot(residuals(NPSA1))

# Respiration
# Eliminate missing R data
RData <- StressMaster[complete.cases(StressMaster$tfRSA),]

# Model selection
RSA1 = lmer(tfRSA ~ Origin * Transplant * Species * Treatment + (1|Parent), data = RData)
RSA2 = lmer(log(tfRSA) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = RData)
RSA3 = lmer(sqrt(tfRSA) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = RData)

model.sel(RSA1, RSA2, RSA3)
summary(RSA1)
anova(RSA1)
shapiro.test(residuals(RSA1))
qqPlot(residuals(RSA1))

# P:R
StressMaster$P_R <- StressMaster$tfGPSA/StressMaster$tfRSA
PRData <- StressMaster[complete.cases(StressMaster$P_R),]

# Model selection
P_R1 = lmer(P_R ~ Origin * Transplant * Species * Treatment + (1|Parent), data = PRData)
P_R2 = lmer(log(P_R) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = PRData)
P_R3 = lmer(sqrt(P_R) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = PRData)

model.sel(P_R1, P_R2, P_R3)
summary(P_R1)
anova(P_R1)
shapiro.test(residuals(P_R3))
qqPlot(residuals(P_R3))

# Calcification
BWData <- StressMaster[complete.cases(StressMaster$mg.cm2.d),]

# Model selection
BWSA1 = lmer(mg.cm2.d ~ Origin * Transplant * Species * Treatment + (1|Parent), data = BWData)
BWSA2 = lmer(log(mg.cm2.d) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = BWData)
BWSA3 = lmer(sqrt(mg.cm2.d) ~ Origin * Transplant * Species * Treatment + (1|Parent), data = BWData)

model.sel(BWSA1, BWSA2, BWSA3)
summary(BWSA1)
anova(BWSA1)
shapiro.test(residuals(BWSA1))
qqPlot(residuals(BWSA1))
