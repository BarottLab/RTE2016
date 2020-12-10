rm(list=ls()) # removes all prior objects

library("reshape") 
library("plotrix")
library("ggplot2")
library("reshape2")
library("dplyr")
library("ggThemeAssist")
library("plyr")
library("pscl")
library("gridExtra")
library("car")
library("lsmeans")
library("nlme")
library("multcomp")
library("multcompView")
library("coefplot")
library("gridExtra")

#####Sedimentation Rate#####
####Load Dataset####
seds<-read.csv("Data/sedimentation.csv", header=TRUE, sep=",", na.strings="NA")
seds$Trap.ID<-as.factor(seds$Trap.ID)
seds$Site.ID<-as.factor(seds$Site.ID)
seds$Rack.ID<-as.factor(seds$Rack.ID)
seds$Time<-as.factor(seds$Time)
seds$Date<-as.factor(seds$Date)
seds$Week<-as.factor(seds$Week)


####Sedimentation by Site####
Seds.Summarize.Site <- ddply(seds, c("Site.ID"), summarise,
                             N    = length(Sediment.Day[!is.na(Sediment.Day)]),
                             mean = mean(Sediment.Day, na.rm=TRUE),
                             sd   = sd(Sediment.Day, na.rm=TRUE),
                             se   = sd / sqrt(N)
);Seds.Summarize.Site

meansite=c(Seds.Summarize.Site$mean)
semsite=c(Seds.Summarize.Site$se)

plotCI(barplot(height=meansite, margin(t=1, r=1, b=1, l=1, unit="pt"), mgp=c(1.7,0.6,0), names.arg=Seds.Summarize.Site$Site.ID, ylim=c(0,0.6), beside = TRUE, main = "Sedimentation Rate", xlab="Site ID", ylab="Sedimentation (g) per day", xpd = FALSE, axes = TRUE), y=meansite, uiw=semsite, liw = semsite, ui=NULL, li=NULL, err="y", add=TRUE, pch=NA, gap=0)

hist(seds$Sediment.Day)
shapiro.test(seds$Sediment.Day)

#non-normal, use non-parametric Wilcoxon Rank Sum Test (by site)

wilcox.test(Sediment.Day~Site.ID, data=seds)
#p<0.001