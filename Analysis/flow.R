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

####Load Dataset####
flow<-read.csv("Data/flow.csv", header=TRUE, sep=",", na.strings="NA")
flow$Site.ID<-as.factor(flow$Site.ID)
flow$Rack.ID<-as.factor(flow$Rack.ID)
flow$Time<-as.factor(flow$Time)
flow$Date<-as.factor(flow$Date)
flow$Week<-as.factor(flow$Week)

####Flow by Site.ID####
Flow.Summarize.Site <- ddply(flow, c("Site.ID"), summarise,
                             N    = length(Percent.Dissolution.Hr[!is.na(Percent.Dissolution.Hr)]),
                             mean = mean(Percent.Dissolution.Hr, na.rm=TRUE),
                             sd   = sd(Percent.Dissolution.Hr, na.rm=TRUE),
                             se   = sd / sqrt(N)
);Flow.Summarize.Site

meansitef=c(Flow.Summarize.Site$mean)
semsitef=c(Flow.Summarize.Site$se)

plotCI(barplot(height=meansitef, margin(t=1, r=1, b=1, l=1, unit="pt"), mgp=c(1.7,0.6,0), names.arg=Flow.Summarize.Site$Site.ID, ylim=c(0,3), beside = TRUE, main = "Relative Water Flow", xlab="Patch Reef", ylab=expression(paste("Dissolution ","(% ","hour"^-1, ")")), xpd = FALSE, axes = TRUE), y=meansitef, uiw=semsitef, liw = semsitef, ui=NULL, li=NULL, err="y", add=TRUE, pch=NA, gap=0, scol="black")

hist(flow$Percent.Dissolution.Hr)
shapiro.test(flow$Percent.Dissolution.Hr)

#non-normal distribution, use a non parametric wilcoson rank sum test (by site)
wilcox.test(Percent.Dissolution.Hr~Site.ID, data=flow)
#significant effect of site (W=16, p<0.001)
