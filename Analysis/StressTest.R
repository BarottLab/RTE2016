library(ggplot2)
library(patchwork)

master <- read.csv("Data/Stress_Master.csv")
master$Transplant <- gsub("Inner Bay", "Inner Lagoon", master$Transplant)
master$Transplant <- gsub("Outer Bay", "Outer Lagoon", master$Transplant)
master$Origin <- gsub("Inner Bay", "Inner Lagoon", master$Origin)
master$Origin <- gsub("Outer Bay", "Outer Lagoon", master$Origin)
master$Species <- gsub("MCAP", "Montipora capitata", as.character(master$Species))
master$Species <- gsub("PCOM", "Porites compressa", as.character(master$Species))
master$PR <- master$GPSA/master$RSA
master$GPSAhour <- master$GPSA/24*1440
master$RSAhour <- master$RSA/24*1440

t0 <- subset(master, Timepoint == 0)
t0am <- subset(t0, Treatment == "Ambient")
t0h <- subset(t0, Treatment == "High")
tf <- subset(master, Timepoint == "f")
tfam <- subset(tf, Treatment == "Ambient")
tfh <- subset(tf, Treatment == "High")

GP0 = ggplot(t0am, aes(x = Transplant, y = GPSAhour, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = t0h, aes(x = Transplant, y = GPSAhour, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = t0h, aes(x = Transplant, y = GPSAhour, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = t0h, aes(x = Transplant, y = GPSAhour, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Gross Photosynthesis", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "", color = "Origin", title = "", tag = "B") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  theme(strip.text.y = element_blank()) +
  coord_cartesian(ylim = c(0, 4)) 
GP0

GPf = ggplot(tfam, aes(x = Transplant, y = GPSAhour, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = tfh, aes(x = Transplant, y = GPSAhour, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = GPSAhour, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = GPSAhour, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Gross Photosynthesis", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "Destination", color = "Origin", title = "", tag = "B") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  theme(strip.text.y = element_blank()) +
  coord_cartesian(ylim = c(0, 4)) 
GPf

R0 = ggplot(t0am, aes(x = Transplant, y = RSAhour, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = t0h, aes(x = Transplant, y = RSAhour, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = t0h, aes(x = Transplant, y = RSAhour, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = t0h, aes(x = Transplant, y = RSAhour, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19), guide = F) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Respiration", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "", color = "Origin", title = "", tag = "C") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  theme(strip.text.y = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5)) 
R0

Rf = ggplot(tfam, aes(x = Transplant, y = RSAhour, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = tfh, aes(x = Transplant, y = RSAhour, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = RSAhour, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = RSAhour, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Respiration", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = "Destination", color = "Origin", title = "", tag = "C") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  theme(strip.text.y = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5)) 
Rf

BWf = ggplot(tfam, aes(x = Transplant, y = x.d, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = tfh, aes(x = Transplant, y = x.d, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = x.d, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = x.d, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Calcification", paste(("%"~day^-1)))), x = "Destination", color = "Origin", title = "", tag = "D") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  theme(strip.text.y = element_blank())
BWf

BWtf = ggplot(tfam, aes(x = Transplant, y = mg.cm2.d, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = tfh, aes(x = Transplant, y = x.d, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = x.d, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = x.d, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Calcification", paste((mg~cm^-2~day^-1)))), x = "Destination", color = "Origin", title = "", tag = "D") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  theme(strip.text.y = element_blank())
BWtf

PAM0 = ggplot(t0am, aes(x = Transplant, y = Yield, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = t0h, aes(x = Transplant, y = Yield, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = t0h, aes(x = Transplant, y = Yield, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = t0h, aes(x = Transplant, y = Yield, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Dark Adapted Yield", paste((F["v"]/F["m"])))), x = "", color = "Origin", title = "", tag = "A") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_text(face = "italic", size = 10), strip.background = element_blank()) +
  theme(strip.text.y = element_blank())
PAM0

PAMf = ggplot(tfam, aes(x = Transplant, y = Yield, color = Origin)) + 
  stat_summary(aes(group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group = Origin), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8, position = position_dodge(width = 0.2)) +  
  stat_summary(data = tfh, aes(x = Transplant, y = Yield, color = Origin, group = Origin, shape = Treatment), fun = mean, geom = "point", size = 2, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = Yield, color = Origin, group = Origin, linetype = Treatment), fun = mean, geom = "line", size = 0.8, position = position_dodge(width = 0.2)) +
  stat_summary(data = tfh, aes(x = Transplant, y = Yield, color = Origin, group = Origin), fun.data = mean_se, geom = "errorbar", size = 0.8, width = 0.2, position = position_dodge(width = 0.2)) +
  scale_x_discrete(expand = c(0.17, 0.17)) +
  scale_shape_manual(values=c("High" = 1, "Ambient" = 19)) +
  scale_linetype_manual(values=c("High" = "dashed", "Ambient" = "solid")) +
  scale_color_manual(values = c("Inner Lagoon" = "coral2", "Outer Lagoon" = "skyblue4")) + 
  theme(aspect.ratio = 1, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10, face = "italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "right", legend.key = element_rect(fill = NA), legend.key.width = unit(0.75,"cm"),
        axis.title.x = element_blank(), plot.tag = element_text(face = "bold"), plot.tag.position = c(0.13, 1)) + 
  scale_fill_manual(values=alpha(c("coral2", "skyblue4")), name = "Origin Site", labels = c("Inner Lagoon", "Outer Lagoon")) + 
  labs(y = expression(atop("Dark Adapted Yield", paste((F["v"]/F["m"])))), x = "Destination", color = "Origin", title = "", tag = "A") + 
  facet_grid(cols = vars(Species)) + 
  theme(strip.text.x = element_text(face = "italic", size = 10), strip.background = element_blank()) +
  theme(strip.text.y = element_blank()) +
  guides(colour = guide_legend(override.aes = list(shape = NA))) 
PAMf



