#Sedimentation Analyses#
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

#####Sedimentation Rate#####
####Load Dataset####
seds<-read.csv("~/Desktop/RTE - PNAS/sedimentation.csv", header=TRUE, sep=",", na.strings="NA")
seds$Trap.ID<-as.factor(seds$Trap.ID)
seds$Site.ID<-as.factor(seds$Site.ID)
seds$Rack.ID<-as.factor(seds$Rack.ID)
seds$Time<-as.factor(seds$Time)
seds$Date<-as.factor(seds$Date)
seds$Week<-as.factor(seds$Week)


####Sedimentation by Site.ID####
Seds.Summarize.Site <- ddply(seds, c("Site.ID"), summarise,
                             N    = length(Sediment.Day),
                             mean = mean(Sediment.Day, na.rm=TRUE),
                             sd   = sd(Sediment.Day, na.rm=TRUE),
                             se   = sd / sqrt(N)
)
Seds.Summarize.Site

meansite=c(Seds.Summarize.Site$mean)
semsite=c(Seds.Summarize.Site$se)

plotCI(barplot(height=meansite, margin(t=1, r=1, b=1, l=1, unit="pt"), mgp=c(1.7,0.6,0), names.arg=Seds.Summarize.Site$Site.ID, ylim=c(0,0.6), beside = TRUE, main = "Sedimentation Rate", xlab="Site ID", ylab="Sedimentation (g) per day", xpd = FALSE, axes = TRUE), y=meansite, uiw=semsite, liw = semsite, ui=NULL, li=NULL, err="y", add=TRUE, pch=NA, gap=0)

hist(seds$Sediment.Day)
shapiro.test(seds$Sediment.Day)
#non-normal, use non-parametric?

wilcox.test(Sediment.Day~Site.ID, alternative=c("greater"), data=seds)

Seds.Site.lm <- lm(Sediment.Day ~ Site.ID, data=seds)
anova(Seds.Site.lm)

Seds.Site.posthoc <- lsmeans(Seds.Site.lm, specs=c("Site.ID"), na.rm=TRUE) #calculate MS means
Seds.Site.posthoc #view results
Seds.Site.posthoc.lett <- cld(Seds.Site.posthoc , alpha=.05, Letters=letters) #identify posthoc letter differences
Seds.Site.posthoc.lett

####Sedimentation by Trap.ID####
Seds.Summarize <- ddply(seds, c("Rack.ID"), summarise,
                        N    = length(Sediment.Day),
                        mean = mean(Sediment.Day, na.rm=TRUE),
                        sd   = sd(Sediment.Day, na.rm=TRUE),
                        se   = sd / sqrt(N)
)
Seds.Summarize

meanstrap=c(Seds.Summarize$mean)
semtrap=c(Seds.Summarize$se)

plotCI(barplot(height=meanstrap, names.arg=Seds.Summarize$Rack.ID, ylim=c(0,1), beside = TRUE, main = "Sedimentation per day by frame", xlab="Frame ID", ylab="Sedimentation (g) per day", xpd = FALSE, axes = TRUE), y=meanstrap, uiw=semtrap, liw = semtrap, ui=NULL, li=NULL, err="y", add=TRUE, pch=NA, gap=0)

Seds.lm <- lm(Sediment.Day ~ Trap.ID, data=seds)
anova(Seds.lm)

Seds.posthoc <- lsmeans(Seds.lm, specs=c("Trap.ID"), na.rm=TRUE) #calculate MS means
Seds.posthoc #view results
Seds.posthoc.lett <- cld(Seds.posthoc , alpha=.05, Letters=letters) #identify posthoc letter differences
Seds.posthoc.lett


####Sedimentation over Time####

Seds.Time<- ddply(seds, c("Week", "Site.ID"), summarise,
                          N    = length(Sediment.Day),
                          mean = mean(Sediment.Day, na.rm=TRUE),
                          sd   = sd(Sediment.Day, na.rm=TRUE),
                          se   = sd / sqrt(N)
)
Seds.Time

Fig1<-ggplot(data=Seds.Time, aes(x=Week, y=mean, group=Site.ID, shape=Site.ID, color=Site.ID))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0, position=position_dodge()) +
  geom_line() +
  geom_point()+
  xlab("Time (weeks)")+
  ylab(expression(Sedimentation~(g~day^{-1})))+
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
        plot.title=element_text(), #Justify the title to the top left
        legend.text = element_text(size = 20),
        axis.title.y=element_text(),
        text=element_text(size=20))
Fig1


#####Sediment Grain Size#####
####Load Dataset####
size<-read.csv("./Sedimentation/SedimentAnalysis.csv", header=TRUE, sep=",", na.strings="NA")
size$Date<-as.POSIXct(size$Date, format="%m/%d/%Y")

#look at overal proportion of sediments in each size class at each site
bins<-ddply(size, c("Site", "Bin"), summarise,
                 N    = length(Proportion),
                 mean = mean(Proportion, na.rm=TRUE),
                 sd   = sd(Proportion, na.rm=TRUE),
                 se   = sd / sqrt(N)
)
bins<-bins[-c(9,10,11,12,13,22,23,24,25,26), ]
bins

binsbar<-ggplot(data=bins, aes(x=Bin, y=mean, fill=Site)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0, position=position_dodge(0.9))+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14))+
  ylab(expression(paste("Proportion of Grains in Size Bin"))) +
  xlab("Size Bins (mm)")
binsbar

#look at mean sediment size at each site

#subset data to pull out only mean values
mean<-size[size$Bin == 'Mean',]
mean

hist(mean$Count)
shapiro.test(mean$Count)
#normally distributed

means<-ddply(mean, c("Site"), summarise,
            N    = length(Count),
            mean = mean(Count, na.rm=TRUE),
            sd   = sd(Count, na.rm=TRUE),
            se   = sd / sqrt(N)
)
means

meansbar<-ggplot(data=means, aes(x=Site, y=mean, fill=Site)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0, position=position_dodge(0.9))+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14))+
  ylab(expression(paste("Mean Size (mm)"))) +
  xlab("Site")
meansbar

summary(aov(mean$Count~mean$Site))
#not significantly different by site

summary(aov(mean$Count~mean$Site*mean$Date))
#when you look at effects of site*date, there is a significant effect of site, date, and interaction of the two


#look at size of mean sediment over time periods (x = date, color = site, line = mean)

meantime<-ddply(mean, c("Date", "Site"), summarise,
                N    = length(Count),
                mean = mean(Count, na.rm=TRUE),
                sd   = sd(Count, na.rm=TRUE),
                se   = sd / sqrt(N)
)
meantime

FigMeanTime<-ggplot(data=meantime, aes(x=Date, y=mean, group=Site, color=Site))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0, position=position_dodge(0)) +
  geom_line() +
  geom_point()+
  xlab("Date")+
  ylab("Mean Size (mm)")+
  theme_bw()+
  theme(axis.title.x=element_text(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_text(),
        legend.position=("right"),
        plot.title=element_text(), #Justify the title to the top left
        legend.text = element_text(size = 20),
        axis.title.y=element_text(),
        text=element_text(size=20))
FigMeanTime




