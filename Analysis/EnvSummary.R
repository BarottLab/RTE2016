## RTE Env Analysis - PNAS ##
setwd("~/Desktop/RTE - PNAS")
library(readr)
library(lubridate)
library(ggplot2)
library(plyr)

# Import data
rm(list = ls())
seaphox <- read_csv("seaphox_final.csv") 
CocoTemp <- read_csv("Coco.csv")
MontyTemp <- read_csv("Monty.csv")

# Format seaphox data
seaphox$newTime <- as.POSIXct(seaphox$Time, format = "%H:%M") # convert time to POSIXct format
seaphox$Time.r <- lubridate::round_date(seaphox$newTime, "15 minutes") # group times by nearest 15-min interval
seaphox$Date <- as.Date(seaphox$Date, format = "%m/%d/%y") # convert date to Date format
seaphox$TimeFinal <- format(seaphox$Time.r, "%H:%M") # convert time format to H:M
seaphox$dttm <- as.POSIXct(paste(seaphox$Date, seaphox$TimeFinal), format="%Y-%m-%d %H:%M") # merge date and time
seaphox$dttm <- as.POSIXlt(seaphox$dttm, format = "%H:%M")
attr(seaphox$dttm, "tzone") <- "UTC"
seaphox$dttm <- as.POSIXct(seaphox$dttm)
seaphox$dttm <- format(seaphox$dttm, tz = "HST")

# Format CTD data
CocoTemp$Site <- "Coco 4"
MontyTemp$Site <- "Monty 13"
Sites <- rbind(CocoTemp, MontyTemp) 
Sites$newTime <- paste(floor(Sites$Time), round((Sites$Time-floor(Sites$Time))*60), sep = ":")
Sites$newTime <- as.POSIXct(Sites$newTime, format = "%H:%M") # convert time to POSIXct format
Sites$Time.r <- lubridate::round_date(Sites$newTime, "15 minutes") # group times by nearest 15-min interval
Sites$newDate <- as.Date(Sites$Day, origin = "2016-01-01") # convert date to Date format
Sites$TimeFinal <- format(Sites$Time.r, "%H:%M") # convert time format to H:M
Sites$dttm <- as.POSIXct(paste(Sites$newDate, Sites$TimeFinal), format="%Y-%m-%d %H:%M") # merge date and time
Sites$dttm <- as.POSIXlt(Sites$dttm, format = "%H:%M")
attr(Sites$dttm, "tzone") <- "UTC"
Sites$dttm <- as.POSIXct(Sites$dttm)
Sites$dttm <- format(Sites$dttm, tz = "HST")

# Merge all data 
Final <- merge(seaphox, Sites, by = "dttm", all = T) # add TempCTD to final seaphox data
Final <- Final[ , c(1:2, 5:8, 19)] # remove unneccessary columns
names(Final) <- c("DTTM", "Site", "pH", "Temp", "Sal", "DO", "TempCTD") # rename columns
Final$DO <- Final$DO/32*1000 # convert DO units
Final <- subset(Final, DTTM < "2017-04-19") # trim by deployment period
Final$DTTM <- as_datetime(Final$DTTM)
Final <- Final[(order(Final$Site)),]

# Separate date and time
Final$DTTM <- as_datetime(Final$DTTM)
Final$Time <- format(Final$DTTM, "%H:%M")
Final$Date <- as.Date(Final$DTTM)
Final$Site <- gsub("Monty 13", "Outer Lagoon", as.character(Final$Site))
Final$Site <- gsub("Coco 4", "Inner Lagoon", as.character(Final$Site))

# Subset rte period
T6 <- subset(Final, DTTM < "2017-02-16" & DTTM > "2016-08-25")
T6 <- T6 %>% distinct(DTTM, Site, .keep_all = TRUE)
T6 <- T6[!is.na(T6$Site),]

#######Load in T6 file - AH ########

T6 <- read.csv("EnvData.csv", na.strings="NA", sep=",")

# Summary Stats
## Temp
Temp.Summarize.Site <- ddply(T6, c("Site"), summarise,
                             N    = length(Temp[!is.na(Temp)]), #AH changed to not include NA's
                             mean = mean(Temp, na.rm=TRUE),
                             sd   = sd(Temp, na.rm=TRUE),
                             se   = sd / sqrt(N)
);Temp.Summarize.Site

ggplot(Temp.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$Temp)
shapiro.test(T6$Temp)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(Temp ~ Site, data = T6) 

# parametric
Temp.Site.lm <- lm(Temp ~ Site, data = T6)
anova(Temp.Site.lm) # p < 0.01

t.test(T6$Temp ~ T6$Site) # p < 0.01

Temp.Site.posthoc <- lsmeans(Temp.Site.lm, specs=c("Site"), na.rm=TRUE) 
Temp.Site.posthoc
Temp.Site.posthoc.lett <- cld(Temp.Site.posthoc , alpha=.05, Letters=letters) 
Temp.Site.posthoc.lett

Temp.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                  N    = length(Temp[!is.na(Temp)]), #AH changed to not include NA's
                                  mean = mean(Temp, na.rm=TRUE)
)

ddply(Temp.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean[!is.na(mean)]), #AH changed to not include NA's
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

Temp.Summarize.Site.Date <- ddply(T6, c("Site", "Date"), summarise,
                                  N    = length(Temp[!is.na(Temp)]), #AH changed to not include NA's
                                  mean = mean(Temp, na.rm=TRUE),
                                  sd   = sd(Temp, na.rm=TRUE)
);Temp.Summarize.Site.Date

wilcox.test(mean ~ Site, data = Temp.Summarize.Site.Date) # p = 0.0948

ggplot() +
  geom_line(subset(Temp.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(Temp.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(Temp.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(Temp.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  theme(panel.background = element_blank()) +
  labs(y = "mean Temp")

## pH
pH.Summarize.Site <- ddply(T6, c("Site"), summarise,
                             N    = length(pH[!is.na(pH)]), #AH changed to not include NA's
                             mean = mean(pH, na.rm=TRUE),
                             sd   = sd(pH, na.rm=TRUE),
                             se   = sd / sqrt(N)
);pH.Summarize.Site

ggplot(pH.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$pH)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(pH ~ Site, data = T6) 

# parametric
pH.Site.lm <- lm(pH ~ Site, data = T6)
anova(pH.Site.lm) # p < 0.01

t.test(T6$pH ~ T6$Site) # p < 0.01

pH.Site.posthoc <- lsmeans(pH.Site.lm, specs=c("Site"), na.rm=TRUE) 
pH.Site.posthoc
pH.Site.posthoc.lett <- cld(pH.Site.posthoc , alpha=.05, Letters=letters) 
pH.Site.posthoc.lett

pH.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                  N    = length(pH[!is.na(pH)]), #AH changed to not include NA's
                                  mean = mean(pH, na.rm=TRUE)
)

ddply(pH.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean[!is.na(mean)]), #AH changed to not include NA's
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

pH.Summarize.Site.Date <- ddply(T6, c("Site", "Date"), summarise,
                                  N    = length(pH[!is.na(pH)]), #AH changed to not include NA's
                                  mean = mean(pH, na.rm=TRUE),
                                sd   = sd(pH, na.rm=TRUE)
);pH.Summarize.Site.Date

wilcox.test(mean ~ Site, data = pH.Summarize.Site.Date) # p < 0.01

ggplot() +
  geom_line(subset(pH.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(pH.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(pH.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(pH.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  theme(panel.background = element_blank()) +
  labs(y = "mean pH")

## Sal
Sal.Summarize.Site <- ddply(T6, c("Site"), summarise,
                           N    = length(Sal[!is.na(Sal)]), #AH changed to not include NA's
                           mean = mean(Sal, na.rm=TRUE),
                           sd   = sd(Sal, na.rm=TRUE),
                           se   = sd / sqrt(N)
);Sal.Summarize.Site

ggplot(Sal.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$Sal)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(Sal ~ Site, data = T6) 

# parametric
Sal.Site.lm <- lm(Sal ~ Site, data = T6)
anova(Sal.Site.lm) # p < 0.01

t.test(T6$Sal ~ T6$Site) # p < 0.01

Sal.Site.posthoc <- lsmeans(Sal.Site.lm, specs=c("Site"), na.rm=TRUE) 
Sal.Site.posthoc
Sal.Site.posthoc.lett <- cld(Sal.Site.posthoc , alpha=.05, Letters=letters) 
Sal.Site.posthoc.lett

Sal.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                N    = length(Sal[!is.na(Sal)]), #AH changed to not include NA's
                                mean = mean(Sal, na.rm=TRUE)
)

ddply(Sal.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean[!is.na(mean)]), #AH changed to not include NA's
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

Sal.Summarize.Site.Date <- ddply(T6, c("Site", "Date"), summarise,
                                N    = length(Sal[!is.na(Sal)]), #AH changed to not include NA's
                                mean = mean(Sal, na.rm=TRUE),
                                sd   = sd(Sal, na.rm=TRUE)
);Sal.Summarize.Site.Date

wilcox.test(mean ~ Site, data = Sal.Summarize.Site.Date) # p < 0.01

ggplot() +
  geom_line(subset(Sal.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(Sal.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(Sal.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(Sal.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  theme(panel.background = element_blank()) +
  labs(y = "mean Sal")

## DO
DO.Summarize.Site <- ddply(T6, c("Site"), summarise,
                            N    = length(DO[!is.na(DO)]), #AH changed to not include NA's,
                            mean = mean(DO, na.rm=TRUE),
                            sd   = sd(DO, na.rm=TRUE),
                            se   = sd / sqrt(N)
);DO.Summarize.Site

ggplot(DO.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$DO)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(DO ~ Site, data = T6) # p < 0.01

# parametric
DO.Site.lm <- lm(DO ~ Site, data = T6)
anova(DO.Site.lm) # p < 0.01

t.test(T6$DO ~ T6$Site) # p < 0.01

DO.Site.posthoc <- lsmeans(DO.Site.lm, specs=c("Site"), na.rm=TRUE) 
DO.Site.posthoc
DO.Site.posthoc.lett <- cld(DO.Site.posthoc , alpha=.05, Letters=letters) 
DO.Site.posthoc.lett

DO.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                N    = length(DO[!is.na(DO)]), #AH changed to not include NA's
                                mean = mean(DO, na.rm=TRUE)
)

ddply(DO.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean[!is.na(mean)]), #AH changed to not include NA's
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

DO.Summarize.Site.Date <- ddply(T6, c("Site", "Date"), summarise,
                                N    = length(DO[!is.na(DO)]), #AH changed to not include NA's
                                mean = mean(DO, na.rm=TRUE)
);DO.Summarize.Site.Date

wilcox.test(mean ~ Site, data = DO.Summarize.Site.Date) # p < 0.01

DO.Summarize.Site.Date <- ddply(T6, c("Site", "Date"), summarise,
                                 N    = length(DO[!is.na(DO)]), #AH changed to not include NA's
                                 mean = mean(DO, na.rm=TRUE),
                                 sd   = sd(DO, na.rm=TRUE)
);DO.Summarize.Site.Date

wilcox.test(mean ~ Site, data = DO.Summarize.Site.Date) # p < 0.01

ggplot() +
  geom_line(subset(DO.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(DO.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(DO.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(DO.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  theme(panel.background = element_blank()) +
  labs(y = "mean DO")

