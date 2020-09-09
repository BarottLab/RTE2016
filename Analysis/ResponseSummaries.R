# Multi-panel Response Figures
## Load necessary packages
library(ggplot2) # plotting
library(patchwork) # arrange plots
library(ggpubr)

## Import dataset
Master <- read.csv("Data/RTE_Frags.csv")

## T6 ##
Tfinal <- subset(Master, Timepoint == "6")
Tfinal$Species <- gsub("MCAP", "Montipora capitata", as.character(Tfinal$Species))
Tfinal$Species <- gsub("PCOM", "Porites compressa", as.character(Tfinal$Species))
Tfinal$Transplant <- gsub("Inner Bay", "Inner Lagoon", Tfinal$Transplant)
Tfinal$Transplant <- gsub("Outer Bay", "Outer Lagoon", Tfinal$Transplant)
Tfinal$Origin <- gsub("Inner Bay", "Inner Lagoon", Tfinal$Origin)
Tfinal$Origin <- gsub("Outer Bay", "Outer Lagoon", Tfinal$Origin)

# Biomass
BMfinal <- subset(Tfinal, BM >= 0)
BM = ggplot(BMfinal, aes(x = Transplant, y = BM, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 10)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y =  element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Biomass", paste((mg~cm^-2)))), x = "Destination", color = "Origin ", title = "", tag = "A") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_text(face = "italic", size = 10), strip.background = element_blank())
BM

# Gross Photosynthesis
Tfinal$GPhour <- Tfinal$GPSA/24*1440

GP = ggplot(Tfinal, aes(x = Transplant, y = GPhour, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 3)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Gross Photosynthesis", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "Destination", color = "Origin ", title = "", tag = "B") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_blank(), strip.background = element_blank())
GP

# Respiration
Tfinal$Rhour <- Tfinal$RSA/24*1440

R = ggplot(Tfinal, aes(x = Transplant, y = Rhour, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Respiration", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "Destination", color = "Origin ", title = "", tag = "C") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_blank(), strip.background = element_blank())
R

# P:R
PR = ggplot(Tfinal, aes(x = Transplant, y = P_RSA, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(2, 4)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title.x=element_text(size=10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 9, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = "P:R", x = "Destination", color = "Origin ", title = "", tag = "D") + 
  facet_grid(cols = vars(Species)) +  theme(strip.text.x = element_blank())
PR

# Survival
SurvivalFinal <- read.csv("RTE_Parents_Master_Long.csv")
SurvivalFinal$Species <- gsub("MCAP", "Montipora capitata", as.character(SurvivalFinal$Species))
SurvivalFinal$Species <- gsub("PCOM", "Porites compressa", as.character(SurvivalFinal$Species))
SurvivalFinal <- subset(SurvivalFinal, Timepoint == "6")
SurvivalFinal$Transplant <- gsub("Inner Bay", "Inner Lagoon", SurvivalFinal$Transplant)
SurvivalFinal$Transplant <- gsub("Outer Bay", "Outer Lagoon", SurvivalFinal$Transplant)
SurvivalFinal$Origin <- gsub("Inner Bay", "Inner Lagoon", SurvivalFinal$Origin)
SurvivalFinal$Origin <- gsub("Outer Bay", "Outer Lagoon", SurvivalFinal$Origin)

S = ggplot(SurvivalFinal, aes(x = Transplant, y = Survival, color = Origin)) + 
  geom_line(aes(group = Parent), size = 0.4, alpha = 0.25) + 
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(50, 100)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 9, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Survival", paste(('%')))), x = "Destination", color = "Origin ", title = "", tag = "E") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_text(face = "italic", size = 10), strip.background = element_blank())
S

# Growth Rate
BWdatafinal <- subset(Tfinal, BWx.d > 0)

BW = ggplot(BWdatafinal, aes(x = Transplant, y = BWx.d, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Calcification", paste(('%'~day^-1)))), x = "Destination", color = "Origin ", title = "", tag = "F") + 
  facet_grid(cols = vars(Species)) +  theme(strip.text.x = element_blank())
BW

# Linear Extension
LEdatafinal <- subset(Tfinal, LE >= 0)
LEdatafinal$LEmm <- LEdatafinal$LE*10
LE = ggplot(LEdatafinal, aes(x = Transplant, y = LEmm, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Linear Extension", paste((mm~day^-1)))), x = "Destination", color = "Origin ", title = "", tag = "G") + 
  facet_grid(cols = vars(Species)) +  theme(strip.text.x = element_blank())
LE

# Reproduction
Repro <- read.csv("RTE_Parents_Master_Long.csv")
Repro$Species <- gsub("MCAP", "Montipora capitata", as.character(Repro$Species))
Repro$Species <- gsub("PCOM", "Porites compressa", as.character(Repro$Species))
Repro <- subset(Repro, Species == "Montipora capitata")
Repro <- subset(Repro, Timepoint == 6)
Repro$Transplant <- gsub("Inner Bay", "Inner Lagoon", Repro$Transplant)
Repro$Transplant <- gsub("Outer Bay", "Outer Lagoon", Repro$Transplant)
Repro$Origin <- gsub("Inner Bay", "Inner Lagoon", Repro$Origin)
Repro$Origin <- gsub("Outer Bay", "Outer Lagoon", Repro$Origin)

Reproduction = ggplot(Repro, aes(x = Transplant, y = PropSpawn, color = Origin)) + 
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(size = 9, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_text(size=10, hjust = 0.5, face = "bold"), legend.key = element_blank(), legend.text = element_text(size = 10), legend.position = c(1.5,0.5), legend.key.size = unit(0.6, "cm"),
        legend.background = element_rect(fill = "white", linetype = "solid", colour = "black", size = 0.25), axis.title.x = element_blank(), 
        plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Reproductive Success", paste((Proportion~Spawned)))), x = "Destination", color = "Origin", title = "", tag = "H")
Reproduction

## Arrange Figures
BM + S + GP + BW + R + LE + PR + {Reproduction + plot_spacer()} + plot_layout(ncol = 2, nrow = 4, heights = 2)

## T3 ##
Tmid <- subset(Master, Timepoint == "3")
Tmid$Species <- gsub("MCAP", "Montipora capitata", as.character(Tmid$Species))
Tmid$Species <- gsub("PCOM", "Porites compressa", as.character(Tmid$Species))
Tmid$Transplant <- gsub("Inner Bay", "Inner Lagoon", Tmid$Transplant)
Tmid$Transplant <- gsub("Outer Bay", "Outer Lagoon", Tmid$Transplant)
Tmid$Origin <- gsub("Inner Bay", "Inner Lagoon", Tmid$Origin)
Tmid$Origin <- gsub("Outer Bay", "Outer Lagoon", Tmid$Origin)

# Biomass
BMdatamid <- subset(Tmid, BM >= 0)
BM = ggplot(BMdatamid, aes(x = Transplant, y = BM, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 15)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y =  element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Biomass", paste((mg~cm^-2)))), x = "Destination", color = "Origin ", title = "", tag = "A") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_text(face = "italic", size = 10), strip.background = element_blank())
BM

# Gross Photosynthesis
Tmid$GPhour <- Tmid$GPSA/24*1440

GP = ggplot(Tmid, aes(x = Transplant, y = GPhour, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 3)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Gross Photosynthesis", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "Destination", color = "Origin ", title = "", tag = "B") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_blank(), strip.background = element_blank())
GP

# Respiration
Tmid$Rhour <- Tmid$RSA/24*1440

R = ggplot(Tmid, aes(x = Transplant, y = Rhour, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Respiration", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "Destination", color = "Origin ", title = "", tag = "C") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_blank(), strip.background = element_blank())
R

# P:R
PRmid <- subset(Tmid, P_RSA < 10)
PR = ggplot(PRmid, aes(x = Transplant, y = P_RSA, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(2, 4)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title.x=element_text(size=10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 9, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = "P:R", x = "Destination", color = "Origin ", title = "", tag = "D") + 
  facet_grid(cols = vars(Species)) +  theme(strip.text.x = element_blank())
PR

# Survival
SurvivalMid <- read.csv("RTE_Parents_Master_Long.csv")
SurvivalMid$Species <- gsub("MCAP", "Montipora capitata", as.character(SurvivalMid$Species))
SurvivalMid$Species <- gsub("PCOM", "Porites compressa", as.character(SurvivalMid$Species))
SurvivalMid <- subset(SurvivalMid, Timepoint == "3")
SurvivalMid$Transplant <- gsub("Inner Bay", "Inner Lagoon", SurvivalMid$Transplant)
SurvivalMid$Transplant <- gsub("Outer Bay", "Outer Lagoon", SurvivalMid$Transplant)
SurvivalMid$Origin <- gsub("Inner Bay", "Inner Lagoon", SurvivalMid$Origin)
SurvivalMid$Origin <- gsub("Outer Bay", "Outer Lagoon", SurvivalMid$Origin)

S = ggplot(SurvivalMid, aes(x = Transplant, y = Survival, color = Origin)) + 
  geom_line(aes(group = Parent), size = 0.4, alpha = 0.25) + 
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(50, 100)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 9, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Survival", paste(('%')))), x = "Destination", color = "Origin ", title = "", tag = "E") + 
  facet_grid(cols = vars(Species)) + theme(strip.text.x = element_text(face = "italic", size = 10), strip.background = element_blank())
S

# Growth Rate
BWdatamid <- subset(Tmid, BWx.d > 0)

BW = ggplot(BWdatamid, aes(x = Transplant, y = BWx.d, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Calcification", paste(('%'~day^-1)))), x = "Destination", color = "Origin ", title = "", tag = "F") + 
  facet_grid(cols = vars(Species)) +  theme(strip.text.x = element_blank())
BW

# Linear Extension
LEdatamid <- subset(Tmid, LE >= 0)
LEdatamid$LEmm <- LEdatamid$LE*10
LE = ggplot(LEdatamid, aes(x = Transplant, y = LEmm, color = Origin)) + 
  stat_summary(aes(group = Parent), fun = mean, geom = "line", size = 0.4, alpha = 0.25) +
  stat_summary(aes(group = Origin), fun = mean, geom = "point", size = 1.25, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) + 
  scale_x_discrete(expand = c(0.17, 0.17)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text.x=element_text(size=9, color = "black"), axis.text.y = element_text(size = 9, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin ", labels = c("Inner Lagoon Origin", "Outer Lagoon Origin")) + 
  labs(y = expression(atop("Linear Extension", paste((mm~day^-1)))), x = "Destination", color = "Origin ", title = "", tag = "G") + 
  facet_grid(cols = vars(Species)) +  theme(strip.text.x = element_blank())
LE

## Arrange Figures
BM + S + GP + BW + R + LE + PR + plot_layout(ncol = 2, nrow = 4, heights = 2)