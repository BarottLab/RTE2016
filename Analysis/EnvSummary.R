setwd("~/Documents/Professional/GitHub/RTE2016")
library(readr)
library(lubridate)
library(ggplot2)
library(plyr)

# Import data
T6 <- read.csv("Data/EnvData.csv", na.strings="NA", sep=",")

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

Temp.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                  N    = length(Temp[!is.na(Temp)]), #AH changed to not include NA's
                                  mean = mean(Temp, na.rm=TRUE),
                                  sd   = sd(Temp, na.rm=TRUE),
                                  se   = sd / sqrt(N)
);Temp.Summarize.Site.Time

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

pH.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                  N    = length(pH[!is.na(pH)]), #AH changed to not include NA's
                                  mean = mean(pH, na.rm=TRUE)
);pH.Summarize.Site.Time

ddply(pH.Summarize.Site.Time, c("Site"), summarise,
      N = length(mean[!is.na(mean)]), #AH changed to not include NA's
      range = max(mean, na.rm = T) - min(mean, na.rm = T),
      sd = sd(mean, na.rm = T))

pH.Summarize.Site.Date <- ddply(T6, c("Site", "Date"), summarise,
                                N = length(pH[!is.na(pH)]), #AH changed to not include NA's
                                mean = mean(pH, na.rm=TRUE),
                                sd = sd(pH, na.rm=TRUE)
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

Sal.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                N    = length(Sal[!is.na(Sal)]), #AH changed to not include NA's
                                mean = mean(Sal, na.rm=TRUE)
);Sal.Summarize.Site.Time

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

DO.Summarize.Site.Time <- ddply(T6, c("Site", "Time"), summarise,
                                N    = length(DO[!is.na(DO)]), #AH changed to not include NA's
                                mean = mean(DO, na.rm=TRUE)
);DO.Summarize.Site.Time

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
t.test(mean ~ Site, data = DO.Summarize.Site.Date) # p < 0.01
DO.Site.lm <- lm(mean ~ Site+Date, data = DO.Summarize.Site.Date)
anova(DO.Site.lm)

ggplot() +
  geom_line(subset(DO.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(DO.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, y = mean, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(DO.Summarize.Site.Date, Site == "Outer Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(DO.Summarize.Site.Date, Site == "Inner Lagoon"), mapping = aes(x = Date, ymin = mean-sd, ymax = mean+sd, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  theme(panel.background = element_blank()) +
  labs(y = "mean DO")

