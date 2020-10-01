#Master RTE Analysis Script for Coral Responses#

####Setup####
rm(list=ls(all=TRUE)) 
library("plotrix") #functions in tapply
library("ggplot2") #plotting
library("reshape2") #reshape data
library("dplyr")  #splitting, applying, and combining data
library("plyr")  #splitting, applying, and combining data
library("pscl")
library("rcompanion")
library("plotrix") #plotting
library("car") #levenes test
library("lsmeans")
library("effects")
library("multcomp") 
library("vegan") #calculating distance matrices
library("betareg")
library("coin")
library("glmmADMB")
library("blmeco")
library("PMCMR")
library("MuMIn")
library("lme4")
library("tidyverse")
library("cowplot")
library("lmerTest")
library("ggfortify")
library("vegan")
library("glmmTMB")
library("emmeans")

####Load Datasets####
#fragment data
frags<-read.csv("Data/RTE_Frags.csv", sep=",", na.strings="NA")
#data summarized by genotype/parent
parents<-read.csv("Data/RTE_Frags.csv", sep=",", na.strings="NA")
#spawning data (M. capitata)
spawning<-read.csv("Data/SpawnDynamic.csv", sep=",", na.strings="NA")

#structure frag dataframe
frags$FragID<-as.factor(frags$FragID)
frags$Parent<-as.factor(frags$Parent)
frags$Timepoint<-as.factor(frags$Timepoint)

#structure parent data frame
parents$Parent<-as.factor(parents$Parent)
parents$Timepoint<-as.factor(parents$Timepoint)
#add in species for parent dataframe
parents$Species<-frags$Species[match(parents$Parent, frags$Parent)]

#structure spawning dataframe
spawning$Night<-as.factor(spawning$Night)
spawning$SiteOrigin<-as.factor(spawning$SiteOrigin)
spawning$SiteTransplant<-as.factor(spawning$SiteTransplant)

#conduct transformations prior to subsetting - see analyses below for comparison of transformations
frags$tBiomass<-log(frags$BM) #biomass log transform
frags$tFeed<-sqrt(7+frags$Final.prey.SA) #add constant for negative values
frags$PropLipids<-frags$PercentLipids/100 #calculate percent lipids as proportion lipids

#subset by timepoint and species
mid <- frags[which(frags$Timepoint == "3"),]
final <- frags[which(frags$Timepoint == "6"),]
mid1.5<- frags[which(frags$Timepoint == "1.5"),]
mid4.5<- frags[which(frags$Timepoint == "4.5"),]
pcom <- frags[which(frags$Species == "PCOM"),]
mcap <- frags[which(frags$Species == "MCAP"),]
pcom3 <- pcom[which(pcom$Timepoint == "3"),]
mcap3 <- mcap[which(mcap$Timepoint == "3"),]
pcom6 <- pcom[which(pcom$Timepoint == "6"),]
mcap6 <- mcap[which(mcap$Timepoint == "6"),]
pcom1.5 <- pcom[which(pcom$Timepoint == "1.5"),]
mcap1.5 <- mcap[which(mcap$Timepoint == "1.5"),]
pcom4.5 <- pcom[which(pcom$Timepoint == "4.5"),]
mcap4.5 <- mcap[which(mcap$Timepoint == "4.5"),]
parentsMCAP <- parents[which(parents$Species == "MCAP"),]
parentsPCOM <- parents[which(parents$Species == "PCOM"),]
parents6 <- parents[which(parents$Timepoint == "3"),]
parents3 <- parents[which(parents$Timepoint == "6"),]

####Biomass####

#create data frame for biomass analysis
biomassmid<-mid[complete.cases(mid$tBiomass), ]
biomassmid$group<-paste(biomassmid$Species, biomassmid$Origin, biomassmid$Transplant)
biomassfinal<-final[complete.cases(final$tBiomass), ]
biomassfinal$group<-paste(biomassfinal$Species, biomassfinal$Origin, biomassfinal$Transplant)

#summarize means
ddply(biomassmid, c("Species", "Origin", "Transplant"), summarise,
             N    = length(BM[!is.na(BM)]),
             mean = mean(BM, na.rm=TRUE),
             sd   = sd(BM, na.rm=TRUE),
             se   = sd / sqrt(N)
)

ddply(biomassfinal, c("Species", "Origin", "Transplant"), summarise,
      N    = length(BM[!is.na(BM)]),
      mean = mean(BM, na.rm=TRUE),
      sd   = sd(BM, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model at 3-months for both species

#check model fit with transformation (tBiomass was log transformmed above)
#log transformed
biomass1 = lmer(tBiomass ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=biomassmid)
qqPlot(residuals(biomass1))

#no transformation
biomass2 = lmer((BM) ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=biomassmid)
qqPlot(residuals(biomass2))

#squareroot transformation
biomass3 = lmer(sqrt(BM) ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=biomassmid)
qqPlot(residuals(biomass3))

model.sel(biomass1, biomass2, biomass3)
#log transformation is best fit according to delta AIC

summary(biomass1)
leveneTest(residuals(biomass1)~biomassmid$group)
anova(biomass1, type=3)

#conduct post hoc test
emm = emmeans(biomass1, ~ Transplant * Origin | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model at 6 months both species 
biomass4 = lmer(tBiomass ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), data=biomassfinal)
qqPlot(residuals(biomass4))
summary(biomass4)
anova(biomass4, type=3)
leveneTest(residuals(biomass4)~biomassfinal$group)

#conduct post hoc test
emm = emmeans(biomass4, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))




####Heterotrophic Feeding#### 

#create data frame for feeding analysis
feedingmid<-mid[complete.cases(mid$tFeed), ]
feedingmid$group<-paste(feedingmid$Species, feedingmid$Origin, feedingmid$Transplant)
feedingfinal<-final[complete.cases(final$tFeed), ]
feedingfinal$group<-paste(feedingfinal$Species, feedingfinal$Origin, feedingfinal$Transplant)

#summarize means
ddply(feedingmid, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Final.prey.SA[!is.na(Final.prey.SA)]),
      mean = mean(Final.prey.SA, na.rm=TRUE),
      sd   = sd(Final.prey.SA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(feedingfinal, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Final.prey.SA[!is.na(Final.prey.SA)]),
      mean = mean(Final.prey.SA, na.rm=TRUE),
      sd   = sd(Final.prey.SA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model 3 months both species
#sqrt transformed
feeding1 = lmer(tFeed ~ Origin * Transplant * Species + (1|Parent), data=feedingmid)
qqPlot(residuals(feeding1))

#no transformation
feeding2 = lmer(Final.prey.SA ~ Origin * Transplant * Species + (1|Parent), data=feedingmid)
qqPlot(residuals(feeding2))

#compare models
model.sel(feeding1, feeding2)
#square root transformed data set has lowest AIC score

summary(feeding1)
leveneTest(residuals(feeding1)~feedingmid$group)
anova(feeding1, type=3)

#conduct post hoc test
emm = emmeans(feeding1, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model 6 months both species
feeding3 = lmer(tFeed ~ Origin * Transplant * Species + (1|Parent), data=feedingfinal)
qqPlot(residuals(feeding3))
summary(feeding3)
anova(feeding3, type=3)
leveneTest(residuals(feeding3)~feedingfinal$group)

#conduct post hoc test
emm = emmeans(feeding3, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

####Lipid Content####

#create data frame for lipid analysis

lipidfinal<-final[complete.cases(final$PropLipids), ]
lipidfinal$group<-paste(lipidfinal$Species, lipidfinal$Origin, lipidfinal$Transplant)

#summarize means
ddply(lipidfinal, c("Species", "Origin", "Transplant"), summarise,
      N    = length(PropLipids[!is.na(PropLipids)]),
      mean = mean(PropLipids, na.rm=TRUE),
      sd   = sd(PropLipids, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#modeling proportion lipids (between 0 and 1) with a logit link function in the glmmTMB package with beta error structure
library(glmmadmb)
lipids1<-glmmTMB(PropLipids ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), lipidfinal, beta_family())
qqPlot(residuals(lipids1))
summary(lipids1)
Anova(lipids1, type="III") 
leveneTest(lipidfinal$PropLipids~lipidfinal$group)

#conduct post hoc test
emm = emmeans(lipids1, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))


####Proportion Spawned####

spawning$ParentID <- ifelse(spawning$ColonyID %in% c("21A", "21B"), 21,
                     ifelse(spawning$ColonyID %in% c("22A", "22B"), 22, 
                     ifelse(spawning$ColonyID %in% c("23A", "23B"), 23,
                     ifelse(spawning$ColonyID %in% c("24A", "24B"), 24,
                     ifelse(spawning$ColonyID %in% c("51A", "51B"), 51,
                     ifelse(spawning$ColonyID %in% c("52A", "52B"), 52,
                     ifelse(spawning$ColonyID%in% c("55A", "55B"), 55,
                     ifelse(spawning$ColonyID %in% c("143A", "143B"), 143,
                     ifelse(spawning$ColonyID %in% c("144A", "144B"), 144,
                     ifelse(spawning$ColonyID %in% c("208A", "208B"), 208,
                     ifelse(spawning$ColonyID %in% c("71A", "71B"), 71,
                     ifelse(spawning$ColonyID %in% c("72A", "72B"), 72,
                     ifelse(spawning$ColonyID %in% c("121A", "121B"), 121,
                     ifelse(spawning$ColonyID %in% c("125A", "125B"), 125,
                     ifelse(spawning$ColonyID %in% c("162A", "162B"), 162,
                     ifelse(spawning$ColonyID %in% c("163A", "163B"), 163,
                     ifelse(spawning$ColonyID %in% c("166A", "166B"), 166,
                     ifelse(spawning$ColonyID %in% c("167A", "167B"), 167,
                     ifelse(spawning$ColonyID %in% c("169A", "169B"), 169,
                     ifelse(spawning$ColonyID %in% c("170A", "170B"), 170, NA))))))))))))))))))))

spawning$Origin <- ifelse(spawning$SiteOrigin %in% 4, "Inner",
                          ifelse(spawning$SiteOrigin %in% 13, "Outer", NA))

spawning$Transplant <- ifelse(spawning$SiteTransplant %in% 4, "Inner",
                              ifelse(spawning$SiteTransplant %in% 13, "Outer", NA))

#summarize means
spawningmeans<-ddply(spawning, c("Origin", "Transplant", "ColonyID"), summarise,
      N    = length(Spawn[!is.na(Spawn)]),
      mean = mean(Spawn, na.rm=TRUE),
      sd   = sd(Spawn, na.rm=TRUE),
      se   = sd / sqrt(N)
);spawningmeans

ddply(spawningmeans, c("Origin", "Transplant"), summarise,
      N    = length(mean[!is.na(mean)]),
      mean = mean(mean, na.rm=TRUE),
      sd   = sd(mean, na.rm=TRUE),
      se   = sd / sqrt(N)
)


#conduct model
spawn.mod <- glmer(Spawn~Origin*Transplant + (1|ParentID) + (1|Date), family = binomial, spawning)

summary(spawn.mod)
Anova(spawn.mod, type = 3)
emmeans(spawn.mod, pairwise~Origin | Transplant, adjust = "tukey")
dispersion_glmer(spawn.mod) #below recommended threshold of 1.4, no evidence of overdispersion

#run overdispersion function from B. Bolker
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(spawn.mod) 


####PAM#### 

pam1.5<-mid1.5[complete.cases(mid1.5$Yield), ]
pam1.5$group<-paste(pam1.5$Species, pam1.5$Origin, pam1.5$Transplant)

pam3<-mid[complete.cases(mid$Yield), ]
pam3$group<-paste(pam3$Species, pam3$Origin, pam3$Transplant)

pam4.5<-mid4.5[complete.cases(mid4.5$Yield), ]
pam4.5$group<-paste(pam4.5$Species, pam4.5$Origin, pam4.5$Transplant)

#summarize means
ddply(pam1.5, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Yield[!is.na(Yield)]),
      mean = mean(Yield, na.rm=TRUE),
      sd   = sd(Yield, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(pam3, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Yield[!is.na(Yield)]),
      mean = mean(Yield, na.rm=TRUE),
      sd   = sd(Yield, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(pam4.5, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Yield[!is.na(Yield)]),
      mean = mean(Yield, na.rm=TRUE),
      sd   = sd(Yield, na.rm=TRUE),
      se   = sd / sqrt(N)
)


#model 1.5 months both species

#test for transformations
pam1a = lmer(Yield ~ Origin * Transplant * Species + (1|Parent)  + (1|Rack), data=pam1.5)
qqPlot(residuals(pam1a))

pam1b = lmer(log(Yield) ~ Origin * Transplant * Species + (1|Parent)  + (1|Rack), data=pam1.5)
qqPlot(residuals(pam1b))

pam1c = lmer((Yield)^2 ~ Origin * Transplant * Species + (1|Parent)  + (1|Rack), data=pam1.5)
qqPlot(residuals(pam1c))

#compare models
model.sel(pam1a, pam1b, pam1c)
#model with log transformation is the best fit according to AIC scores

summary(pam1b)
leveneTest(residuals(pam1b)~pam1.5$group) 
anova(pam1b, type=3)

#conduct post hoc test
emm = emmeans(pam1b, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))


#model 3 months both species
pam2a = lmer(Yield ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=pam3)
qqPlot(residuals(pam2a))

pam2b = lmer(log(Yield) ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=pam3)
qqPlot(residuals(pam2b))

pam2c = lmer((Yield)^2 ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=pam3)
qqPlot(residuals(pam2c))

#compare models
model.sel(pam2a, pam2b, pam2c)
#model with log transformation is the best fit according to AIC scores

summary(pam2b)
leveneTest(residuals(pam2b)~pam3$group)
anova(pam2b, type=3)

#conduct post hoc test
emm = emmeans(pam2b, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))



#model 4.5 months both species
pam3a = lmer(Yield ~ Origin * Transplant * Species + (1|Parent)+ (1|Rack), data=pam4.5)
qqPlot(residuals(pam3a))

pam3b = lmer(log(Yield) ~ Origin * Transplant * Species + (1|Parent)+ (1|Rack), data=pam4.5)
qqPlot(residuals(pam3b))

pam3c = lmer(Yield^2 ~ Origin * Transplant * Species + (1|Parent)+ (1|Rack), data=pam4.5)
qqPlot(residuals(pam3c))

#compare models
model.sel(pam3a, pam3b, pam3c)
#model with log transformation is the best fit according to AIC scores

summary(pam3b)
leveneTest(residuals(pam3b)~pam4.5$group)
anova(pam3b, type=3)

#conduct post hoc test
emm = emmeans(pam3b, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))



####Survival####

survfinal<-parents6[complete.cases(parents6$Survival), ]
survfinal$group<-paste(survfinal$Species, survfinal$Origin, survfinal$Transplant)

survmid<-parents3[complete.cases(parents3$Survival), ]
survmid$group<-paste(survmid$Species, survmid$Origin, survmid$Transplant)


#convert to proportion rather than percent as currently organized
survfinal$Survival<-survfinal$Survival/100
survmid$Survival<-survmid$Survival/100

#summarize means
ddply(survmid, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Survival[!is.na(Survival)]),
      mean = mean(Survival, na.rm=TRUE),
      sd   = sd(Survival, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(survfinal, c("Species", "Origin", "Transplant"), summarise,
      N    = length(Survival[!is.na(Survival)]),
      mean = mean(Survival, na.rm=TRUE),
      sd   = sd(Survival, na.rm=TRUE),
      se   = sd / sqrt(N)
)


#change survival to be 0.9999 for 1's
survfinal$Beta<-ifelse(survfinal$Survival==1, 0.99999, survfinal$Survival)
survmid$Beta<-ifelse(survmid$Survival==1, 0.99999, survmid$Survival)

#analyze 3-month survivorship with a beta regression
surv1<-betareg(Beta ~ Origin * Transplant * Species, data=survmid)
qqPlot(residuals(surv1))
summary(surv1)
Anova(surv1, type="III") #chi square tests

#conduct post hoc test
emm = emmeans(surv1, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#analyze 6-month survivorship
surv2<-betareg(Beta ~ Origin * Transplant * Species, data=survfinal)
qqPlot(residuals(surv2))
summary(surv2)
Anova(surv2, type="III") 

#conduct post hoc test
emm = emmeans(surv2, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

####Fitness Scores####

## with no reproduction 

#create data frame for fitness analysis
GSfitnessMCAP<-parentsMCAP[complete.cases(parentsMCAP$FitScoreNoRepro), ]
GSfitnessMCAP$group<-paste(GSfitnessMCAP$Species, GSfitnessMCAP$Origin, GSfitnessMCAP$Transplant)

GSfitnessPCOM<-parentsPCOM[complete.cases(parentsPCOM$FitScoreNoRepro), ]
GSfitnessPCOM$group<-paste(GSfitnessPCOM$Species, GSfitnessPCOM$Origin, GSfitnessPCOM$Transplant)

#summarize means
ddply(GSfitnessMCAP, c("Species", "Origin", "Transplant"), summarise,
      N    = length(FitScoreNoRepro[!is.na(FitScoreNoRepro)]),
      mean = mean(FitScoreNoRepro, na.rm=TRUE),
      sd   = sd(FitScoreNoRepro, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(GSfitnessPCOM, c("Species", "Origin", "Transplant"), summarise,
      N    = length(FitScoreNoRepro[!is.na(FitScoreNoRepro)]),
      mean = mean(FitScoreNoRepro, na.rm=TRUE),
      sd   = sd(FitScoreNoRepro, na.rm=TRUE),
      se   = sd / sqrt(N)
)

##Model fitness scores for Montipora capitata

fitnessA = lmer(FitScoreNoRepro ~ Origin * Transplant + (1|Parent), data=GSfitnessMCAP)
qqPlot(residuals(fitnessA))
summary(fitnessA)
anova(fitnessA, type=3)
leveneTest(GSfitnessMCAP$FitScoreNoRepro~GSfitnessMCAP$group)
#passes Levene Test for variance

#conduct post hoc test
emm = emmeans(fitnessA, ~ Origin * Transplant, adjust="tukey")
#letter display
cld(emm, Letters=c(LETTERS))

##Model fitness scores for Porites compressa

fitnessB = lmer(FitScoreNoRepro ~ Origin * Transplant + (1|Parent), data=GSfitnessPCOM)
qqPlot(residuals(fitnessB))
summary(fitnessB)
anova(fitnessB, type=3)
leveneTest(GSfitnessPCOM$FitScoreNoRepro~GSfitnessPCOM$group)
#passes Levene Test for variance

#conduct post hoc test
emm = emmeans(fitnessB, ~ Origin * Transplant, adjust="tukey")
#letter display
cld(emm, Letters=c(LETTERS))



## with reproduction

#create data frame for fitness analysis
fitnessMCAP<-parentsMCAP[complete.cases(parentsMCAP$FitScore), ]
fitnessMCAP$group<-paste(fitnessMCAP$Species, fitnessMCAP$Origin, fitnessMCAP$Transplant)

fitnessPCOM<-parentsPCOM[complete.cases(parentsPCOM$FitScore), ]
fitnessPCOM$group<-paste(fitnessPCOM$Species, fitnessPCOM$Origin, fitnessPCOM$Transplant)

#summarize means
ddply(fitnessMCAP, c("Species", "Origin", "Transplant"), summarise,
      N    = length(FitScore[!is.na(FitScore)]),
      mean = mean(FitScore, na.rm=TRUE),
      sd   = sd(FitScore, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(fitnessPCOM, c("Species", "Origin", "Transplant"), summarise,
      N    = length(FitScore[!is.na(FitScore)]),
      mean = mean(FitScore, na.rm=TRUE),
      sd   = sd(FitScore, na.rm=TRUE),
      se   = sd / sqrt(N)
)

##Model fitness scores for Montipora capitata

fitness1 = lmer(FitScore ~ Origin * Transplant + (1|Parent), data=fitnessMCAP)
qqPlot(residuals(fitness1))
summary(fitness1)
anova(fitness1, type=3)
leveneTest(fitnessPCOM$FitScore~fitnessPCOM$group)
#passes Levene Test for variance

#conduct post hoc test
emm = emmeans(fitness1, ~ Origin * Transplant, adjust="tukey")
#letter display
cld(emm, Letters=c(LETTERS))

##Model fitness scores for Porites compressa

fitness2 = lmer(FitScore ~ Origin * Transplant + (1|Parent), data=fitnessPCOM)
qqPlot(residuals(fitness2))
summary(fitness2)
anova(fitness2, type=3)
leveneTest(fitnessPCOM$FitScore~fitnessPCOM$group)
#passes Levene Test for variance

#conduct post hoc test
emm = emmeans(fitness2, ~ Origin * Transplant, adjust="tukey")
#letter display
cld(emm, Letters=c(LETTERS))



####Linear Extension####

#linear extension from long form data frame (use "LE")

#create data frame for linear extension
le03<-mid[complete.cases(mid$LE), ]
le03$group<-paste(le03$Species, le03$Origin, le03$Transplant)
le03 <- subset(le03, LE > 0)
le03$LE <- le03$LE*10

le06<-final[complete.cases(final$LE), ]
le06$group<-paste(le06$Species, le06$Origin, le06$Transplant)
le06 <- subset(le06, LE > 0)
le06$LE <- le06$LE*10

#summarize means
ddply(le03, c("Species", "Origin", "Transplant"), summarise,
      N    = length(LE[!is.na(LE)]),
      mean = mean(LE, na.rm=TRUE),
      sd   = sd(LE, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(le06, c("Species", "Origin", "Transplant"), summarise,
      N    = length(LE[!is.na(LE)]),
      mean = mean(LE, na.rm=TRUE),
      sd   = sd(LE, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model 0-3 months both species
le1 = lmer(sqrt(LE) ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=le03)
qqPlot(residuals(le1))

le1b = lmer(LE ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), data=le03)
qqPlot(residuals(le1b))

#model selection
model.sel(le1, le1b)

#original data is the best fit

summary(le1)
leveneTest(residuals(le1)~le03$group)
anova(le1, type=3) 

#conduct post hoc test
emm = emmeans(le1, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model 0-6 months both species
le2 = lmer(LE ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), data=le06)
qqPlot(residuals(le2))

le2b = lmer(sqrt(LE) ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), data=le06)
qqPlot(residuals(le2b))

#model selection
model.sel(le2, le2b)
#untransformed data is the best

summary(le2)
leveneTest(residuals(le2)~le06$group) 
anova(le2, type=3)

#conduct post hoc test
emm = emmeans(le2, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))


####Buoyant Weight####

#from long form data frame at 3 and 6 month timepoints 

#create data frame for buoyant weight
bw03<-mid[complete.cases(mid$BWx.d), ]
bw03$group<-paste(bw03$Species, bw03$Origin, bw03$Transplant)
bw03 <- subset(bw03, BWx.d > 0)
bw06<-final[complete.cases(final$BWx.d), ]
bw06$group<-paste(bw06$Species, bw06$Origin, bw06$Transplant)
bw06 <- subset(bw06, BWx.d > 0)

#summarize means
ddply(bw03, c("Species", "Origin", "Transplant"), summarise,
      N    = length(BWx.d[!is.na(BWx.d)]),
      mean = mean(BWx.d, na.rm=TRUE),
      sd   = sd(BWx.d, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(bw06, c("Species", "Origin", "Transplant"), summarise,
      N    = length(BWx.d[!is.na(BWx.d)]),
      mean = mean(BWx.d, na.rm=TRUE),
      sd   = sd(BWx.d, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model 0-3 months both species
bw1 = lmer(sqrt(BWx.d) ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), data=bw03) 
#sqrt transformation meets assumptions

summary(bw1)
qqPlot(residuals(bw1))
leveneTest(residuals(bw1)~bw03$group) 
anova(bw1, type=3)

#conduct post hoc test
emm = emmeans(bw1, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))


#model 0-6 months both species
bw2 = lmer(sqrt(BWx.d) ~ Origin * Transplant * Species + (1|Parent) +  (1|Rack), data=bw06)
summary(bw2)
qqPlot(residuals(bw2))
leveneTest(residuals(bw2)~bw06$group)
anova(bw2, type=3)

#conduct post hoc test
emm = emmeans(bw2, ~ Origin * Transplant | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))


####Net Photosynthesis####

#create data frame for net photosynthesis
NP3<-mid[complete.cases(mid$NPSA), ]
NP3$group<-paste(NP3$Species, NP3$Origin, NP3$Transplant)
NP6<-final[complete.cases(final$NPSA), ]
NP6$group<-paste(NP6$Species, NP6$Origin, NP6$Transplant)

#summarize means
ddply(NP3, c("Species", "Origin", "Transplant"), summarise,
      N    = length(NPSA[!is.na(NPSA)]),
      mean = mean(NPSA, na.rm=TRUE),
      sd   = sd(NPSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(NP6, c("Species", "Origin", "Transplant"), summarise,
      N    = length(NPSA[!is.na(NPSA)]),
      mean = mean(NPSA, na.rm=TRUE),
      sd   = sd(NPSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)


#model at 3-months for both species
NPSA1 = lmer(NPSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=NP3)
qqPlot(residuals(NPSA1)) 

summary(NPSA1)
leveneTest(residuals(NPSA1)~NP3$group) 
anova(NPSA1, type=3)

#conduct post hoc test
emm = emmeans(NPSA1, ~ Transplant * Origin * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model at 6 months both species 
NPSA2 = lmer(NPSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=NP6)

qqPlot(residuals(NPSA2)) 
summary(NPSA2)
leveneTest(residuals(NPSA2)~NP6$group) 
anova(NPSA2, type=3)

#conduct post hoc test
emm = emmeans(NPSA2, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))



####Gross Photosynthesis####

#create data frame for gross photosynthesis
GP3<-mid[complete.cases(mid$GPSA), ]
GP3$group<-paste(GP3$Species, GP3$Origin, GP3$Transplant)
GP6<-final[complete.cases(final$GPSA), ]
GP6$group<-paste(GP6$Species, GP6$Origin, GP6$Transplant)

#summarize means
ddply(GP3, c("Species", "Origin", "Transplant"), summarise,
      N    = length(GPSA[!is.na(GPSA)]),
      mean = mean(GPSA, na.rm=TRUE),
      sd   = sd(GPSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(GP6, c("Species", "Origin", "Transplant"), summarise,
      N    = length(GPSA[!is.na(GPSA)]),
      mean = mean(GPSA, na.rm=TRUE),
      sd   = sd(GPSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model at 3-months for both species
GPSA1 = lmer(GPSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=GP3)
qqPlot(residuals(GPSA1)) 
summary(GPSA1)
leveneTest(residuals(GPSA1)~GP3$group) 
anova(GPSA1, type=3)

#conduct post hoc test
emm = emmeans(GPSA1, ~ Transplant * Origin | Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model at 6 months both species 
GPSA2 = lmer(GPSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=GP6)
qqPlot(residuals(GPSA2)) 
summary(GPSA2)
leveneTest(residuals(GPSA2)~GP6$group) 
anova(GPSA2, type=3)

#conduct post hoc test
emm = emmeans(GPSA2, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

####Respiration####

#create data frame for respiration
RS3<-mid[complete.cases(mid$RSA), ]
RS3$group<-paste(RS3$Species, RS3$Origin, RS3$Transplant)
RS6<-final[complete.cases(final$RSA), ]
RS6$group<-paste(RS6$Species, RS6$Origin, RS6$Transplant)

#summarize means
ddply(RS3, c("Species", "Origin", "Transplant"), summarise,
      N    = length(RSA[!is.na(RSA)]),
      mean = mean(RSA, na.rm=TRUE),
      sd   = sd(RSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(RS6, c("Species", "Origin", "Transplant"), summarise,
      N    = length(RSA[!is.na(RSA)]),
      mean = mean(RSA, na.rm=TRUE),
      sd   = sd(RSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model at 3-months for both species
RSA1 = lmer(RSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=RS3)
qqPlot(residuals(RSA1)) 
summary(RSA1)
leveneTest(residuals(RSA1)~RS3$group) 
anova(RSA1, type=3)

#conduct post hoc test
emm = emmeans(RSA1, ~ Transplant * Origin * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model at 6 months both species 
RSA2 = lmer(RSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=RS6)
qqPlot(residuals(RSA2)) 
summary(RSA2)
leveneTest(residuals(RSA2)~RS6$group) 
anova(RSA2, type=3)

#conduct post hoc test
emm = emmeans(RSA2, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

####P:R Ratio#### 

#create data frame for ratio 
PR3<-mid[complete.cases(mid$P_RSA), ]
PR3$group<-paste(PR3$Species, PR3$Origin, PR3$Transplant)
PR6<-final[complete.cases(final$P_RSA), ]
PR6$group<-paste(PR6$Species, PR6$Origin, PR6$Transplant)


#summarize means
ddply(PR3, c("Species", "Origin", "Transplant"), summarise,
      N    = length(P_RSA[!is.na(P_RSA)]),
      mean = mean(P_RSA, na.rm=TRUE),
      sd   = sd(P_RSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

ddply(PR6, c("Species", "Origin", "Transplant"), summarise,
      N    = length(P_RSA[!is.na(P_RSA)]),
      mean = mean(P_RSA, na.rm=TRUE),
      sd   = sd(P_RSA, na.rm=TRUE),
      se   = sd / sqrt(N)
)

#model at 3-months for both species
PR1 = lmer(P_RSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=PR3)
qqPlot(residuals(PR1)) 
summary(PR1)
leveneTest(residuals(PR1)~PR3$group) 
anova(PR1, type=3)

#conduct post hoc test
emm = emmeans(PR1, ~ Transplant * Origin * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))

#model at 6 months both species 
PR2 = lmer(P_RSA ~ Origin * Transplant * Species + (1|Parent) + (1|Rack), data=PR6)
qqPlot(residuals(PR2)) 
summary(PR2)
leveneTest(residuals(PR2)~PR6$group) 
anova(PR2, type=3)

#conduct post hoc test
emm = emmeans(PR2, ~ Origin * Transplant * Species, adjust="tukey")
cld(emm, Letters=c(LETTERS))



