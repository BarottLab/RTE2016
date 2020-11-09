#Clod Card Flow Analyses#
#Last modified 15November2016 by Ariana Huffmyer#
#RTE Field Physical Data Analysis#

####Load workspace####
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

setwd("~/Google Drive/Coral AE Project/Field/FieldPhysicalData/")

####Load Dataset####
flow<-read.csv("flow.csv", header=TRUE, sep=",", na.strings="NA")
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


####Sedimentation by Rack.ID####
Flow.Summarize.Rack <- ddply(flow, c("Rack.ID"), summarise,
                        N    = length(Percent.Dissolution.Hr[!is.na(Percent.Dissolution.Hr)]),
                        mean = mean(Percent.Dissolution.Hr, na.rm=TRUE),
                        sd   = sd(Percent.Dissolution.Hr, na.rm=TRUE),
                        se   = sd / sqrt(N)
);Flow.Summarize.Rack

meanrackf=c(Flow.Summarize.Rack$mean)
semrackf=c(Flow.Summarize.Rack$se)

plotCI(barplot(height=meanrackf, names.arg=Flow.Summarize.Rack$Rack.ID, ylim=c(0,3), beside = TRUE, main = "Dissolution by frame", xlab="Frame ID", ylab="Dissolution % per hour", xpd = FALSE, axes = TRUE), y=meanrackf, uiw=semrackf, liw = semrackf, ui=NULL, li=NULL, err="y", add=TRUE, pch=NA, gap=0)

#similar flow by rack at each site. 

####Sedimentation over Time####

Flow.Time<- ddply(flow, c("Week", "Site.ID"), summarise,
                  N    = length(Percent.Dissolution.Hr[!is.na(Percent.Dissolution.Hr)]),
                  mean = mean(Percent.Dissolution.Hr, na.rm=TRUE),
                  sd   = sd(Percent.Dissolution.Hr, na.rm=TRUE),
                  se   = sd / sqrt(N)
);Flow.Time

Fig1<-ggplot(data=Flow.Time, aes(x=Week, y=mean, group=Site.ID, shape=Site.ID, color=Site.ID))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0, position=position_dodge()) +
  geom_line() +
  geom_point()+
  xlab("Time (Weeks)")+
  ylab(expression(Flow~("%"~dissolution~day^{-1})))+
  ggtitle("")+
  theme_bw()+
  labs(color = "Reef Site", shape="Reef Site")+
  scale_color_manual(values=c("gray", "black"))+
  theme(axis.title.x=element_text(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_text(),
        legend.position=("right"),
        plot.title=element_text(), 
        legend.text = element_text(size = 20),
        axis.title.y=element_text(),
        text=element_text(size=20, color="black"));Fig1

#consistent higher flow in OB site