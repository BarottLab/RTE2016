# Set working directory and load necessary packages
library(reshape2)
library(plyr)
library(ggplot2)
library(lubridate)
library(zoo)
library(plotrix)
library(devtools)
library(tools)
library(readr)
library(dplyr)
library(mgcv)
library(ggpubr)
library(ggrepel)

# Import data and subset by timepoint
PAR <- read.csv("Data/PAR.csv")
RTE <- subset(PAR, Date < "2017-02-17")
SP <- subset(PAR, Date > "2017-02-16")

MeanRTE <- cbind(aggregate(RTE$cal, by = list(RTE$Time, RTE$Site), FUN = mean), aggregate(RTE$cal, by = list(RTE$Time, RTE$Site), FUN = sem), aggregate(RTE$cal, by = list(RTE$Time, RTE$Site), FUN = ci))
MeanRTE <- MeanRTE[,c(1:3,6,9)]
names(MeanRTE) <- c("Time", "Site", "PAR", "PARSEM", "PARCI")

MeanSP <- cbind(aggregate(SP$cal, by = list(SP$Time, SP$Site), FUN = mean), aggregate(SP$cal, by = list(SP$Time, SP$Site), FUN = sem), aggregate(SP$cal, by = list(SP$Time, SP$Site), FUN = ci))
MeanSP <- MeanSP[,c(1:3,6,9)]
names(MeanSP) <- c("Time", "Site", "PAR", "PARSEM", "PARCI")

aggregate(MeanRTE$PAR, by = MeanRTE["Site"], FUN = max) # IB = 528.012, OB = 579.977
aggregate(MeanRTE$PARCI, by = MeanRTE["Site"], FUN = mean) # IB = 8.859, OB = 9.444

# Plot total average PAR per site
T06 <- c("T[0-6]")
T710 <- c("T[7-10]")

AvePAR = ggplot() +
  geom_line(subset(MeanRTE, Site == "OB"), mapping = aes(x = Time, y = PAR, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(MeanRTE, Site == "IB"), mapping = aes(x = Time, y = PAR, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(MeanRTE, Site == "OB"), mapping = aes(x = Time, ymin = PAR-PARCI, ymax = PAR+PARCI, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(MeanRTE, Site == "IB"), mapping = aes(x = Time, ymin = PAR-PARCI, ymax = PAR+PARCI, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  labs(y = "PAR", x = "Time") +
  scale_x_discrete(breaks = c("06:00", "12:00", "18:00")) +
  theme(aspect.ratio = 0.8, axis.text=element_text(size=10), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") 
AvePAR

SPPAR = ggplot() +
  geom_line(subset(MeanSP, Site == "OB"), mapping = aes(x = Time, y = PAR, group = 1), color = "skyblue4", size = 0.5) +
  geom_line(subset(MeanSP, Site == "IB"), mapping = aes(x = Time, y = PAR, group = 1), color = "coral2", size = 0.5) +
  geom_ribbon(subset(MeanSP, Site == "OB"), mapping = aes(x = Time, ymin = PAR-PARCI, ymax = PAR+PARCI, group = 1), fill = "skyblue4", alpha = 0.2, linetype = 4) +
  geom_ribbon(subset(MeanSP, Site == "IB"), mapping = aes(x = Time, ymin = PAR-PARCI, ymax = PAR+PARCI, group = 1), fill = "coral2", alpha = 0.2, linetype = 4) +
  labs(y = "PAR", x = "Time") +
  scale_x_discrete(breaks = c("06:00", "12:00", "18:00")) +
  theme(aspect.ratio = 0.8, axis.text=element_text(size=10), axis.title=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") 
SPPAR

# Calculate daily light integral (DLI)
dli <- rbind(F1, F2, F3, F4, F5, F6)
dli$Site <- ifelse(dli$Rack == "F1" | dli$Rack == "F2" | dli$Rack == "F3", "IB", "OB")
dli$DLI <- dli$cal*0.0864
dlit6 <- subset(dli, Date < "2017-02-17")
meandli <- cbind(aggregate(dli$DLI, by = list(dli$Date, dli$Site), FUN = mean), aggregate(dli$DLI, by = list(dli$Date, dli$Site), FUN = sem), aggregate(dli$DLI, by = list(dli$Date, dli$Site), FUN = ci))
meandli <- meandli[,c(1:3,6,9)]
names(meandli) <- c("Date", "Site", "DLI", "DLISEM", "DLICI")

idx <- c(1, diff(dli$Date))
i2 <- c(1,which(idx != 1), nrow(dli)+1)
dli$grp <- rep(1:length(diff(i2)), diff(i2))

idx <- c(1, diff(meandli$Date))
i2 <- c(1,which(idx != 1), nrow(meandli)+1)
meandli$grp <- rep(1:length(diff(i2)), diff(i2))

# Plot mean DLI
T1 <- c("T[0]")
T3 <- c("T[3]")
T6 <- c("T[6]")
S1 <- c("S[1]")
S2 <- c("S[2]")
S3 <- c("S[3]")

TimeseriesDLI = ggplot(meandli, aes(x = Date, y = DLI, color = Site, fill = Site, group = grp)) + 
  geom_rect(aes(xmin = as_date("2017-05-24"), xmax = as_date("2017-05-28"), ymin = -Inf, ymax = Inf), fill = "grey80", color = "grey80", alpha = 0.1) +
  geom_rect(aes(xmin = as_date("2017-06-23"), xmax = as_date("2017-06-27"), ymin = -Inf, ymax = Inf), fill = "grey80", color = "grey80", alpha = 0.1) +
  geom_rect(aes(xmin = as_date("2017-07-22"), xmax = as_date("2017-07-26"), ymin = -Inf, ymax = Inf), fill = "grey80", color = "grey80", alpha = 0.1) +
  geom_line(size = 0.35) +
  scale_x_date(date_labels = "%b '%y", date_breaks = "2 month", limits = as.Date(c("2016-08-01", "2017-08-10"))) + 
  labs(x = "Date", y = expression(atop("DLI", paste((mol~m^-2~day^-1))))) +
  scale_y_continuous(limits = c(0,40)) +
  geom_vline(xintercept = as_date("2016-08-25"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_date("2016-11-28"), linetype = 2, color = "black") +
  geom_vline(xintercept = as_date("2017-02-16"), linetype = 2, color = "black") +
  annotate("text", x = as_date("2016-09-02"), y = 38, label = T1, colour = "black", size = 3, fontface = 1, parse = TRUE) +
  annotate("text", x = as_date("2016-12-05"), y = 38, label = T3, colour = "black", size = 3, fontface = 1, parse = TRUE) +
  annotate("text", x = as_date("2017-02-24"), y = 38, label = T6, colour = "black", size = 3, fontface = 1, parse = TRUE) +
  annotate("text", x = as_date("2017-06-06"), y = 38, label = S1, colour = "black", size = 3, fontface = 1, parse = TRUE) +
  annotate("text", x = as_date("2017-07-06"), y = 38, label = S2, colour = "black", size = 3, fontface = 1, parse = TRUE) +
  annotate("text", x = as_date("2017-08-04"), y = 38, label = S3, colour = "black", size = 3, fontface = 1, parse = TRUE) + 
  theme(aspect.ratio = .8, axis.text=element_text(size=8), axis.title=element_text(size=10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, linetype = 1), legend.title = element_blank(), legend.key = element_blank(), legend.text = element_text(size = 10), 
        legend.position = c(0.1395, 0.0725), legend.key.size = unit(3, "mm"), legend.background = element_rect(fill = "white", linetype = "solid", colour = "black", size = 0.25)) + 
  scale_color_manual(values = c("coral2", "skyblue4"))
TimeseriesDLI

# Final Light Figure
Odyssey <- ggarrange(TimeseriesDLI, ggarrange(AvePAR, SPPAR, ncol = 1, nrow = 2, labels = c("B", "C")), ncol = 2, labels = "A", widths = 2:1) 
Odyssey

aggregate(dlit6$DLI, by = dlit6["Site"], FUN = mean, na.rm = T) # IB = 11.802, OB = 13.265
aggregate(dlit6$DLI, by = dlit6["Site"], FUN = sd, na.rm = T) # IB = 19.338, OB = 20.611

