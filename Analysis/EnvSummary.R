## RTE Env Analysis - PNAS ##
setwd("~/Box/Barott lab/Data/RTE2016/RTE2016 seaphox")
library(readr)
library(lubridate)
library(ggplot2)
library(plyr)
library(tidyverse)

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

write.csv(T6, "T6seaphox.csv")

# Summary Stats
## Temp
Temp.Summarize.Site <- ddply(T6, c("Site"), summarise,
                             N    = length(TempCTD),
                             mean = mean(TempCTD, na.rm=TRUE),
                             sd   = sd(TempCTD, na.rm=TRUE),
                             se   = sd / sqrt(N)
)
Temp.Summarize.Site

ggplot(Temp.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$TempCTD)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(TempCTD ~ Site, alternative = c("greater"), data = T6) # p = 1

# parametric
Temp.Site.lm <- lm(Temp ~ Site, data = T6)
anova(Temp.Site.lm) # p < 0.01

t.test(T6$TempCTD ~ T6$Site) # p < 0.01

Temp.Site.posthoc <- lsmeans(Temp.Site.lm, specs=c("Site"), na.rm=TRUE) 
Temp.Site.posthoc
Temp.Site.posthoc.lett <- cld(Temp.Site.posthoc , alpha=.05, Letters=letters) 
Temp.Site.posthoc.lett

Temp.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                             N    = length(TempCTD),
                             mean = mean(TempCTD, na.rm=TRUE)
)

ddply(Temp.Summarize.Site.Time, c("Site"), summarise,
        N = length(mean),
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

## pH
pH.Summarize.Site <- ddply(T6, c("Site"), summarise,
                             N    = length(pH),
                             mean = mean(pH, na.rm=TRUE),
                             sd   = sd(pH, na.rm=TRUE),
                             se   = sd / sqrt(N)
)
pH.Summarize.Site

ggplot(pH.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$pH)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(pH ~ Site, alternative = c("greater"), data = T6) # p < 0.01

# parametric
pH.Site.lm <- lm(pH ~ Site, data = T6)
anova(pH.Site.lm) # p < 0.01

t.test(T6$pH ~ T6$Site) # p < 0.01

pH.Site.posthoc <- lsmeans(pH.Site.lm, specs=c("Site"), na.rm=TRUE) 
pH.Site.posthoc
pH.Site.posthoc.lett <- cld(pH.Site.posthoc , alpha=.05, Letters=letters) 
pH.Site.posthoc.lett

pH.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                  N    = length(pH),
                                  mean = mean(pH, na.rm=TRUE)
)

ddply(pH.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean),
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

## Sal
Sal.Summarize.Site <- ddply(T6, c("Site"), summarise,
                           N    = length(Sal),
                           mean = mean(Sal, na.rm=TRUE),
                           sd   = sd(Sal, na.rm=TRUE),
                           se   = sd / sqrt(N)
)
Sal.Summarize.Site

ggplot(Sal.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$Sal)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(Sal ~ Site, alternative = c("greater"), data = T6) # p = 1

# parametric
Sal.Site.lm <- lm(Sal ~ Site, data = T6)
anova(Sal.Site.lm) # p < 0.01

t.test(T6$Sal ~ T6$Site) # p < 0.01

Sal.Site.posthoc <- lsmeans(Sal.Site.lm, specs=c("Site"), na.rm=TRUE) 
Sal.Site.posthoc
Sal.Site.posthoc.lett <- cld(Sal.Site.posthoc , alpha=.05, Letters=letters) 
Sal.Site.posthoc.lett

Sal.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                N    = length(Sal),
                                mean = mean(Sal, na.rm=TRUE)
)

ddply(Sal.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean),
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

## DO
DO.Summarize.Site <- ddply(T6, c("Site"), summarise,
                            N    = length(DO),
                            mean = mean(DO, na.rm=TRUE),
                            sd   = sd(DO, na.rm=TRUE),
                            se   = sd / sqrt(N)
)
DO.Summarize.Site

ggplot(DO.Summarize.Site) +
  geom_bar(aes(x = Site, y = mean), stat="identity") +
  geom_errorbar(aes(x = Site, ymin = mean - se, ymax = mean + se), width=0.4)

hist(T6$DO)

# non-parametric, but we have lots of observations (>10,000/site)
wilcox.test(DO ~ Site, alternative = c("greater"), data = T6) # p < 0.01

# parametric
DO.Site.lm <- lm(DO ~ Site, data = T6)
anova(DO.Site.lm) # p < 0.01

t.test(T6$DO ~ T6$Site) # p < 0.01

DO.Site.posthoc <- lsmeans(DO.Site.lm, specs=c("Site"), na.rm=TRUE) 
DO.Site.posthoc
DO.Site.posthoc.lett <- cld(DO.Site.posthoc , alpha=.05, Letters=letters) 
DO.Site.posthoc.lett

DO.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                N    = length(DO),
                                mean = mean(DO, na.rm=TRUE)
)

ddply(DO.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean),
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))
