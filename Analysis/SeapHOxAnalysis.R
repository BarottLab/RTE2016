# Load necessary packages
library(readr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

# Import data
rm(list = ls())
seaphox <- read_csv("Data/seaphox_final.csv") 
CocoTemp <- read_csv("Data/Coco.csv")
MontyTemp <- read_csv("Data/Monty.csv")

### Format seaphox data
# Convert time to POSIXct format
seaphox$newTime <- as.POSIXct(seaphox$Time, format = "%H:%M")

# Group times by nearest 15-min interval
seaphox$Time.r <- lubridate::round_date(seaphox$newTime, "15 minutes") 

# Convert date to Date format
seaphox$Date <- as.Date(seaphox$Date, format = "%m/%d/%y")

# Convert time format to H:M
seaphox$TimeFinal <- format(seaphox$Time.r, "%H:%M")

# Merge date and time
seaphox$dttm <- as.POSIXct(paste(seaphox$Date, seaphox$TimeFinal), format="%Y-%m-%d %H:%M")

# Change timezone
seaphox$dttm <- as.POSIXlt(seaphox$dttm, format = "%H:%M")
attr(seaphox$dttm, "tzone") <- "UTC"
seaphox$dttm <- as.POSIXct(seaphox$dttm)
seaphox$dttm <- format(seaphox$dttm, tz = "HST")

### Format CTD Temp Data 
# Bind site data
CocoTemp$Site <- "Coco 4"
MontyTemp$Site <- "Monty 13"
Sites <- rbind(CocoTemp, MontyTemp)

# Convert time to POSIXct format
Sites$newTime <- paste(floor(Sites$Time), round((Sites$Time-floor(Sites$Time))*60), sep = ":")
Sites$newTime <- as.POSIXct(Sites$newTime, format = "%H:%M")

# Group times by nearest 15-min interval
Sites$Time.r <- lubridate::round_date(Sites$newTime, "15 minutes") 

# Convert date to Date format
Sites$newDate <- as.Date(Sites$Day, origin = "2016-01-01")

# Convert time format to H:M
Sites$TimeFinal <- format(Sites$Time.r, "%H:%M")

# Merge date and time
Sites$dttm <- as.POSIXct(paste(Sites$newDate, Sites$TimeFinal), format="%Y-%m-%d %H:%M")

# Change timezone
Sites$dttm <- as.POSIXlt(Sites$dttm, format = "%H:%M")
attr(Sites$dttm, "tzone") <- "UTC"
Sites$dttm <- as.POSIXct(Sites$dttm)
Sites$dttm <- format(Sites$dttm, tz = "HST")

### Merge all data 
# Add TempCTD to final seaphox data
Final <- merge(seaphox, Sites, by = "dttm", all = T)

# Remove unneccessary columns
Final <- Final[ , c(1:2, 5:8, 19)]

# Rename columns
names(Final) <- c("DTTM", "Site", "pH", "Temp", "Sal", "DO", "TempCTD")

# Convert DO units
Final$DO <- Final$DO/32*1000

# Trim by deployment period
Final <- subset(Final, DTTM < "2017-04-19")

# Sort
Final$DTTM <- as_datetime(Final$DTTM)
Final <- Final[(order(Final$Site)),]

# Assign group to continuous data
TimeDiff <- difftime(Final$DTTM, lag(Final$DTTM, default = Final$DTTM[1]), units = "days")
Final$grp <- cumsum(ifelse(TimeDiff>2,1,0))

### Average Day Calculations ###
# Make new time column
Final$DTTM <- as_datetime(Final$DTTM)
Final$Time <- format(Final$DTTM, "%H:%M")

# Make new date column 
Final$Date <- as.Date(Final$DTTM)
Final$Site <- gsub("Monty 13", "Outer Lagoon", as.character(Final$Site))
Final$Site <- gsub("Coco 4", "Inner Lagoon", as.character(Final$Site))

# Subset by timepoint
T3 <- subset(Final, DTTM < "2016-12-02" & DTTM > "2016-08-25")
T6 <- subset(Final, DTTM < "2017-02-16" & DTTM > "2016-08-25")

# Aggregate and merge averages, SEM, SD and 95% CI
sem <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))
ci <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))*qt(0.975, length(x)-1)
range <- function(x) max(x, na.rm = T)-min(x, na.rm = T)

HourlyMean <- cbind(aggregate(T6$pH, by = list(T6$Site, T6$Time), FUN = mean, na.rm = T), aggregate(T6$pH, by = list(T6$Site, T6$Time), FUN = sem), aggregate(T6$pH, by = list(T6$Site, T6$Time), FUN = ci), aggregate(T6$pH, by = list(T6$Site, T6$Time), FUN = sd, na.rm = T),
              aggregate(T6$TempCTD, by = list(T6$Site, T6$Time), FUN = mean, na.rm = T), aggregate(T6$TempCTD, by = list(T6$Site, T6$Time), FUN = sem), aggregate(T6$TempCTD, by = list(T6$Site, T6$Time), FUN = ci),aggregate(T6$TempCTD, by = list(T6$Site, T6$Time), FUN = sd, na.rm = T),
              aggregate(T6$Sal, by = list(T6$Site, T6$Time), FUN = mean, na.rm = T), aggregate(T6$Sal, by = list(T6$Site, T6$Time), FUN = sem), aggregate(T6$Sal, by = list(T6$Site, T6$Time), FUN = ci), aggregate(T6$Sal, by = list(T6$Site, T6$Time), FUN = sd, na.rm = T),
              aggregate(T6$DO, by = list(T6$Site, T6$Time), FUN = mean, na.rm = T), aggregate(T6$DO, by = list(T6$Site, T6$Time), FUN = sem), aggregate(T6$DO, by = list(T6$Site, T6$Time), FUN = ci), aggregate(T6$DO, by = list(T6$Site, T6$Time), FUN = sd, na.rm = T))
HourlyMean <- HourlyMean[, c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48)]
names(HourlyMean) <- c("Site", "Time", "pH", "pHSEM", "pHCI", "pHSD", "Temp", "TempSEM", "TempCI", "TempSD", "Sal", "SalSEM", "SalCI", "SalSD", "DO", "DOSEM", "DOCI", "DOSD")

### Timeseries Plots ###
timepoints = data.frame(date = as_datetime(c("2016-08-27" ,"2016-11-30", "2017-02-18")), timepoint = c("T[0]", "T[3]", "T[6]"))

pH = ggplot() +
  geom_line(subset(Final, Site == "Outer Lagoon"), mapping = aes(x = DTTM, y = pH, group = grp), color = "skyblue4", size = 0.35) +
  geom_line(subset(Final, Site == "Inner Lagoon"), mapping = aes(x = DTTM, y = pH, group = grp), color = "coral2", size = 0.35) +
  scale_x_datetime(labels = date_format("%b '%y"), date_breaks = "1 month", limits = as_datetime(c("2016-07-01","2017-04-19"))) +
  labs(y = "pH", x = "") +
  coord_cartesian(ylim = c(7.6, 8.2)) +
  geom_vline(xintercept = as_datetime("2016-08-25"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2016-11-28"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2017-02-16"), linetype = 2, color = "black") +
  geom_text(data = timepoints, mapping = aes(x = date, y = c(8.133, 8.133, 8.133), label = timepoint), size = 4, vjust = -0.4, hjust = 0, inherit.aes = FALSE, color = "black", parse = TRUE) +
  theme(aspect.ratio = .3, axis.text=element_text(size=10), axis.title=element_text(size=12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.position = "none")
pH

Sal = ggplot() +
  geom_line(subset(Final, Site == "Outer Lagoon"), mapping = aes(x = DTTM, y = Sal, group = grp), color = "skyblue4", size = 0.35) +
  geom_line(subset(Final, Site == "Inner Lagoon"), mapping = aes(x = DTTM, y = Sal, group = grp), color = "coral2", size = 0.35) +
  scale_x_datetime(labels = date_format("%b '%y"), date_breaks = "1 month", limits = as_datetime(c("2016-07-01","2017-04-19"))) +
  labs(y = "Salinity (psu)", x = "") +
  coord_cartesian(ylim = c(27, 36)) +
  geom_vline(xintercept = as_datetime("2016-08-25"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2016-11-28"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2017-02-16"), linetype = 2, color = "black") +
  geom_text(data = timepoints, mapping = aes(x = date, y = c(8.133, 8.133, 8.133), label = timepoint), size = 4, vjust = -0.4, hjust = 0, inherit.aes = FALSE, color = "black", parse = TRUE) +
  theme(aspect.ratio = .3, axis.text=element_text(size=10), axis.title=element_text(size=12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.position = "none")
Sal

DO = ggplot() +
  geom_line(subset(Final, Site == "Outer Lagoon"), mapping = aes(x = DTTM, y = DO, group = grp), color = "skyblue4", size = 0.35) +
  geom_line(subset(Final, Site == "Inner Lagoon"), mapping = aes(x = DTTM, y = DO, group = grp), color = "coral2", size = 0.35) +
  scale_x_datetime(labels = date_format("%b '%y"), date_breaks = "1 month", limits = as_datetime(c("2016-07-01","2017-04-19"))) +
  labs(y = "Dissolved Oxygen (µM)", x = "Date") +
  coord_cartesian(ylim = c(0, 350)) +
  geom_vline(xintercept = as_datetime("2016-08-25"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2016-11-28"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2017-02-16"), linetype = 2, color = "black") +
  geom_text(data = timepoints, mapping = aes(x = date, y = c(8.133, 8.133, 8.133), label = timepoint), size = 4, vjust = -0.4, hjust = 0, inherit.aes = FALSE, color = "black", parse = TRUE) +
  theme(aspect.ratio = .3, axis.text=element_text(size=10), axis.title=element_text(size=12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.position = "none")
DO

Temp = ggplot() +
  geom_line(subset(Final, Site == "Outer Lagoon"), mapping = aes(x = DTTM, y = Temp, group = grp), color = "skyblue4", size = 0.35) +
  geom_line(subset(Final, Site == "Inner Lagoon"), mapping = aes(x = DTTM, y = Temp, group = grp), color = "coral2", size = 0.35) +
  scale_x_datetime(labels = date_format("%b '%y"), date_breaks = "1 month", limits = as_datetime(c("2016-07-01","2017-04-19"))) +
  labs(y = "Temperature (°C)", x = "") +
  coord_cartesian(ylim = c(21, 30)) +
  geom_vline(xintercept = as_datetime("2016-08-25"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2016-11-28"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_datetime("2017-02-16"), linetype = 2, color = "black") +
  geom_text(data = timepoints, mapping = aes(x = date, y = c(8.133, 8.133, 8.133), label = timepoint), size = 4, vjust = -0.4, hjust = 0, inherit.aes = FALSE, color = "black", parse = TRUE) +
  theme(aspect.ratio = .3, axis.text=element_text(size=10), axis.title=element_text(size=12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.position = "none")
Temp

### Daily average plots ###
AvepH = ggplot() +
  geom_line(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, y = pH, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, y = pH, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, ymin = pH-pHCI, ymax = pH+pHCI, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, ymin = pH-pHCI, ymax = pH+pHCI, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  labs(y = "HourlyMean pH", x = "") +
  coord_cartesian(ylim = c(7.92, 8.06)) +
  scale_x_discrete(breaks = c("06:00", "12:00", "18:00")) +
  theme(aspect.ratio = 0.8, axis.text=element_text(size=10), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") 
AvepH

AveTemp = ggplot() +
  geom_line(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, y = Temp, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, y = Temp, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, ymin = Temp-TempCI, ymax = Temp+TempCI, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, ymin = Temp-TempCI, ymax = Temp+TempCI, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  labs(y = "HourlyMean Temperature (°C)", x = "") +
  scale_x_discrete(breaks = c("06:00", "12:00", "18:00")) +
  coord_cartesian(ylim = c(24.5, 26)) +
  theme(aspect.ratio = 0.8, axis.text=element_text(size=10), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") 
AveTemp

AveSal = ggplot() +
  geom_line(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, y = Sal, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, y = Sal, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, ymin = Sal-SalCI, ymax = Sal+SalCI, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, ymin = Sal-SalCI, ymax = Sal+SalCI, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  labs(y = "HourlyMean Salinity (psu)", x = "") +
  scale_x_discrete(breaks = c("06:00", "12:00", "18:00")) +
  scale_y_continuous(breaks = seq(33, 35, by = 1), limits = c(33, 35)) +
  theme(aspect.ratio = 0.8, axis.text=element_text(size=10), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") 
AveSal

AveDO = ggplot() +
  geom_line(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, y = DO, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, y = DO, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(HourlyMean, Site == "Outer Lagoon"), mapping = aes(x = Time, ymin = DO-DOCI, ymax = DO+DOCI, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(HourlyMean, Site == "Inner Lagoon"), mapping = aes(x = Time, ymin = DO-DOCI, ymax = DO+DOCI, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  labs(y = "HourlyMean Dissolved Oxygen (µM)", x = "") +
  scale_x_discrete(breaks = c("06:00", "12:00", "18:00")) +
  theme(aspect.ratio = 0.8, axis.text=element_text(size=10), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") 
AveDO

### Organize plots ###
SeapHOxPanels = ggarrange(Temp, AveTemp, pH, AvepH, Sal, AveSal, DO, AveDO, ncol = 2, nrow = 4, widths = 2:1, heights = 1, align = "v", labels = "AUTO", label.x = c(0.1, 0.25, 0.1, 0.25, 0.1, 0.25, 0.1, 0.25), label.y = 0.97)
SeapHOxPanels

### Summary Stats ###
SiteMeans <- cbind(aggregate(T6$pH, by = list(T6$Site), FUN = mean, na.rm = T), aggregate(T6$pH, by = list(T6$Site), FUN = sd, na.rm = T),
                   aggregate(T6$TempCTD, by = list(T6$Site), FUN = mean, na.rm = T), aggregate(T6$TempCTD, by = list(T6$Site), FUN = sd, na.rm = T),
                   aggregate(T6$Sal, by = list(T6$Site), FUN = mean, na.rm = T), aggregate(T6$Sal, by = list(T6$Site), FUN = sd, na.rm = T),
                   aggregate(T6$DO, by = list(T6$Site), FUN = mean, na.rm = T), aggregate(T6$DO, by = list(T6$Site), FUN = sd, na.rm = T))
SiteMeans <- SiteMeans[ , c(1:2, 4, 6, 8, 10, 12, 14, 16)]
names(SiteMeans) <- c("Site", "pH", "pHSD", "Temp", "TempSD", "Sal", "SalSD", "DO", "DOSD")

SiteRanges <- cbind(aggregate(T6$pH, by = list(T6$Site, T6$Date), FUN = range), aggregate(T6$pH, by = list(T6$Site, T6$Date), FUN = sd, na.rm = T),
                   aggregate(T6$Temp, by = list(T6$Site, T6$Date), FUN = range), aggregate(T6$Temp, by = list(T6$Site, T6$Date), FUN = sd, na.rm = T),
                   aggregate(T6$Sal, by = list(T6$Site, T6$Date), FUN = range), aggregate(T6$Sal, by = list(T6$Site, T6$Date), FUN = sd, na.rm = T),
                   aggregate(T6$DO, by = list(T6$Site, T6$Date), FUN = range), aggregate(T6$DO, by = list(T6$Site, T6$Date), FUN = sd, na.rm = T))
SiteRanges <- SiteRanges[ , c(1:3,6,9,12,15,18,21,24)]
names(SiteRanges) <- c("Site", "Date", "pH", "pHSD", "Temp", "TempSD", "Sal", "SalSD", "DO", "DOSD")

MeanRanges <- cbind(aggregate(SiteRanges$pH, by = list(SiteRanges$Site), FUN = mean, na.rm = T), aggregate(SiteRanges$pHSD, by = list(SiteRanges$Site), FUN = mean, na.rm = T),
                    aggregate(SiteRanges$Temp, by = list(SiteRanges$Site), FUN = mean, na.rm = T), aggregate(SiteRanges$TempSD, by = list(SiteRanges$Site), FUN = mean, na.rm = T),
                    aggregate(SiteRanges$Sal, by = list(SiteRanges$Site), FUN = mean, na.rm = T), aggregate(SiteRanges$SalSD, by = list(SiteRanges$Site), FUN = mean, na.rm = T),
                    aggregate(SiteRanges$DO, by = list(SiteRanges$Site), FUN = mean, na.rm = T), aggregate(SiteRanges$DOSD, by = list(SiteRanges$Site), FUN = mean, na.rm = T))
MeanRanges <- MeanRanges[,c(1:2,4,6,8,10,12,14,16)]
names(MeanRanges) <- c("Site", "pH", "pHSD", "Temp", "TempSD", "Sal", "SalSD", "DO", "DOSD")

## Significance tests
# Means
var.test(T6$pH[T6$Site == "Inner Lagoon"], T6$pH[T6$Site == "Outer Lagoon"]) 
t.test(T6$pH[T6$Site == "Inner Lagoon"], T6$pH[T6$Site == "Outer Lagoon"], var.equal = F) 
var.test(T6$TempCTD[T6$Site == "Inner Lagoon"], T6$TempCTD[T6$Site == "Outer Lagoon"]) 
t.test(T6$TempCTD[T6$Site == "Inner Lagoon"], T6$TempCTD[T6$Site == "Outer Lagoon"], var.equal = T) 
var.test(T6$Sal[T6$Site == "Inner Lagoon"], T6$Sal[T6$Site == "Outer Lagoon"]) 
t.test(T6$Sal[T6$Site == "Inner Lagoon"], T6$Sal[T6$Site == "Outer Lagoon"], var.equal = F) 
var.test(T6$DO[T6$Site == "Inner Lagoon"], T6$DO[T6$Site == "Outer Lagoon"]) 
t.test(T6$DO[T6$Site == "Inner Lagoon"], T6$DO[T6$Site == "Outer Lagoon"], var.equal = F) 

# Amplitudes
var.test(SiteRanges$pH[SiteRanges$Site == "Inner Lagoon"], SiteRanges$pH[SiteRanges$Site == "Outer Lagoon"]) 
t.test(SiteRanges$pH[SiteRanges$Site == "Inner Lagoon"], SiteRanges$pH[SiteRanges$Site == "Outer Lagoon"], var.equal = F) 
var.test(SiteRanges$Temp[SiteRanges$Site == "Inner Lagoon"], SiteRanges$Temp[SiteRanges$Site == "Outer Lagoon"]) 
t.test(SiteRanges$Temp[SiteRanges$Site == "Inner Lagoon"], SiteRanges$Temp[SiteRanges$Site == "Outer Lagoon"], var.equal = T) 
var.test(SiteRanges$Sal[SiteRanges$Site == "Inner Lagoon"], SiteRanges$Sal[SiteRanges$Site == "Outer Lagoon"]) 
t.test(SiteRanges$Sal[SiteRanges$Site == "Inner Lagoon"], SiteRanges$Sal[SiteRanges$Site == "Outer Lagoon"], var.equal = F) 
var.test(SiteRanges$DO[SiteRanges$Site == "Inner Lagoon"], SiteRanges$DO[SiteRanges$Site == "Outer Lagoon"]) 
t.test(SiteRanges$DO[SiteRanges$Site == "Inner Lagoon"], SiteRanges$DO[SiteRanges$Site == "Outer Lagoon"], var.equal = F) 
