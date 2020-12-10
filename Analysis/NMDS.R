# set working directory and load necessary packages
library(vegan)
library(tidyverse)
library(reshape2)

# clear environment
rm(list = ls())

# set seed
set.seed(54321)

# import data
genoslong <- read.csv("Data/RTE_Genotypes.csv") #read in file
genoslong$spe.tra <- paste(genoslong$Species, genoslong$Transplant) #add group name for Species and Transplant Location
genoslong$spe.ori <- paste(genoslong$Species, genoslong$Origin) #add group name for Species and Origin Location
genoslong$group <- paste(genoslong$Species, genoslong$Origin, genoslong$Transplant) #add group name for Species and Origin Location

# subset by timepoint and for response variables for multivariate analysis
t3 <- na.omit(subset(genoslong, Timepoint == "3")[,c("Parent","Transplant","Origin","Timepoint","Species","spe.tra", "spe.ori", "group", "BM","BWx.d","LE","GPSA","RSA","P_RSA","Survival")])
t3$identifier <- rownames(t3)
t6 <- na.omit(subset(genoslong, Timepoint == "6")[,c("Parent","Transplant","Origin","Timepoint","Species","spe.tra", "spe.ori", "group", "BM","BWx.d","LE","GPSA","RSA","P_RSA","Survival")])
t6$identifier <- rownames(t6)
t6 <- t6[-69,]

#check values for outliers
#t3
colMeans(t3[9:15])
par(mfrow=c(2,4))
boxplot(t3$BM)
boxplot(t3$BWx.d)
boxplot(t3$LE)
boxplot(t3$GPSA)
boxplot(t3$RSA)
boxplot(t3$P_RSA)
boxplot(t3$Survival)
dev.off()

#t6
colMeans(t6[9:15])
par(mfrow=c(2,4))
boxplot(t6$BM)
boxplot(t6$BWx.d)
boxplot(t6$LE)
boxplot(t6$GPSA)
boxplot(t6$RSA)
boxplot(t6$P_RSA)
boxplot(t6$Survival)
dev.off()

#Scale and center datasets
t3scaled <- scale(t3[9:15], center = T, scale = T) # scaled variables
t6scaled <- scale(t6[9:15], center = T, scale = T) # scaled variables

#Identify Factors 
fac.t3 <- t3[1:8]
fac.t6 <- t6[1:8]


### Time 3
#PCA
T3.pca.out <- prcomp(t3scaled, center=FALSE, scale=FALSE)
summary(T3.pca.out)
biplot(T3.pca.out)
PC1 <- T3.pca.out$x[,1]
PC2 <- T3.pca.out$x[,2]

# PERMANOVA
## t3
T3.mod <- adonis2(t3scaled ~ Transplant * Origin * Species, data = t3, method = "euclidian") # PERMANOVA
T3.mod
vect3 <- envfit(T3.pca.out, t3[9:15], perm = 1000) #fit physiological vectors onto ordination
vect3df <- as.data.frame(vect3$vectors$arrows * sqrt(vect3$vectors$r))
vect3df$variable <- rownames(vect3df)

#testing for differences in betadispersion by 8 groups of species, origin, and transplant
bdisp.8group <- betadisper(vegdist(t3scaled, method = "euclidian"), fac.t3$group)
anova(bdisp.8group)

#testing for differences in betadispersion by 4 groups of species and transplant
bdisp.4group <- betadisper(vegdist(t3scaled, method = "euclidian"), fac.t3$spe.tra)
anova(bdisp.4group)

#gathering sample info for plotting
pca.df.t3 <- data.frame(T3.pca.out$x[,1], y = T3.pca.out$x[,2],
                         Parent = as.factor(fac.t3$Parent),
                         Origin = as.factor(fac.t3$Origin),
                         Transplant = as.factor(fac.t3$Transplant),
                         Species = as.factor(fac.t3$Species),
                         spe.tra = as.factor(fac.t3$spe.tra),
                         spe.ori = as.factor(fac.t3$spe.ori),
                         group = as.factor(fac.t3$group))

#set colors in this order MCAP Inner Lagoon MCAP Outer Lagoon PCOM Inner Lagoon PCOM Outer Lagoon
groupcolors <- c("cyan", "lightgreen", "green", "darkgreen", 
                 "pink", "purple", "orchid1", "darkorchid4") #set 8 group colors
PCA.colors <- c("lightgreen", "darkgreen", "orchid1", "darkorchid4") #set 4 Origin or Transplant colors
PCA.symbols <- c(16,17) #set two species shapes 16 = Mcap, 17 = Pcomp

#ordering levels
pca.df.t3$origins <- factor(pca.df.t3$spe.ori, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
pca.df.t3$transplants <- factor(pca.df.t3$spe.tra, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
pca.df.t3$group <- factor(pca.df.t3$group, levels = c("MCAP Inner Lagoon Inner Lagoon", "MCAP Inner Lagoon Outer Lagoon", "MCAP Outer Lagoon Outer Lagoon","MCAP Outer Lagoon Inner Lagoon", 
                                                        "PCOM Inner Lagoon Inner Lagoon", "PCOM Inner Lagoon Outer Lagoon", "PCOM Outer Lagoon Outer Lagoon", "PCOM Outer Lagoon Inner Lagoon"))
#setting column names
colnames(pca.df.t3)[1:2] <- c("PC1", "PC2")

pdf("RAnalysis/Output/Time3.PCA.pdf", width=8, height=4)
par(mfrow=c(1,2))
# All Group Ordination
T3.PCA.A.plot <- ordiplot(T3.pca.out, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5), xlab=c("PC1 (35%)"), ylab=c("PC2 (23%)"))
points(T3.PCA.A.plot, "sites", col = groupcolors[pca.df.t3$group], cex=0.8, pch=PCA.symbols[pca.df.t3$Species])
ordihull(T3.PCA.A.plot, groups = pca.df.t3$group, draw = "polygon", alpha = 100, col = groupcolors, border = groupcolors)
par.new = T
plot(vect3, col = "black")
title(main = "A) Time 3, All Groups ")

#Transplant Ordination
T3.PCA.T.plot <- ordiplot(T3.pca.out, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5),xlab=c("PC1 (35%)"), ylab=c("PC2 (23%)"))
points(T3.PCA.T.plot, "sites", col = groupcolors[pca.df.t3$group], cex=0.8, pch=PCA.symbols[pca.df.t3$Species])
ordihull(T3.PCA.T.plot, groups = pca.df.t3$spe.tra, draw = "polygon", alpha = 100, col = PCA.colors, border = PCA.colors)
par.new = T
plot(vect3, col = "black")
title(main = "B) Time 3, Transplant Groups")
legend("topright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
                                 "PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
       pch = c(16,16,16,16,17,17,17,17),
       col = groupcolors, cex=0.4, box.lwd = 0,box.col = "white",bg = "white")
dev.off()

# Time 6
#PCA
T6.pca.out <- prcomp(t6scaled, center=FALSE, scale=FALSE)
summary(T6.pca.out)
biplot(T6.pca.out)
PC1 <- T6.pca.out$x[,1]
PC2 <- T6.pca.out$x[,2]

# PERMANOVA
## t6
T6.mod <- adonis2(t6scaled ~ Transplant * Origin * Species, data = t6, method = "euclidian") # PERMANOVA
T6.mod
vect6 <- envfit(T6.pca.out, t6[9:15], perm = 1000) #fit physiological vectors onto ordination
vect6df <- as.data.frame(vect6$vectors$arrows * sqrt(vect6$vectors$r))
vect6df$variable <- rownames(vect6df)

#PCA
T6.pca.out <- prcomp(t6scaled, center=FALSE, scale=FALSE)
summary(T6.pca.out)
biplot(T6.pca.out)
PC1 <- T6.pca.out$x[,1]
PC2 <- T6.pca.out$x[,2]

#testing for differences in betadispersion by 8 groups of species, origin, and transplant
bdisp.8group <- betadisper(vegdist(t6scaled, method = "euclidian"), fac.t6$group)
anova(bdisp.8group)

#testing for differences in betadispersion by 4 groups of species and transplant
bdisp.4group <- betadisper(vegdist(t6scaled, method = "euclidian"), fac.t6$spe.tra)
anova(bdisp.4group)

#gathering sample info for plotting
pca.df.t6 <- data.frame(T6.pca.out$x[,1], y = T6.pca.out$x[,2],
                        Parent = as.factor(fac.t6$Parent),
                        Origin = as.factor(fac.t6$Origin),
                        Transplant = as.factor(fac.t6$Transplant),
                        Species = as.factor(fac.t6$Species),
                        spe.tra = as.factor(fac.t6$spe.tra),
                        spe.ori = as.factor(fac.t6$spe.ori),
                        group = as.factor(fac.t6$group))

#ordering levels
pca.df.t6$origins <- factor(pca.df.t6$spe.ori, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
pca.df.t6$transplants <- factor(pca.df.t6$spe.tra, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
pca.df.t6$group <- factor(pca.df.t6$group, levels = c("MCAP Inner Lagoon Inner Lagoon", "MCAP Inner Lagoon Outer Lagoon", "MCAP Outer Lagoon Outer Lagoon","MCAP Outer Lagoon Inner Lagoon", 
                                                      "PCOM Inner Lagoon Inner Lagoon", "PCOM Inner Lagoon Outer Lagoon", "PCOM Outer Lagoon Outer Lagoon", "PCOM Outer Lagoon Inner Lagoon"))
#setting column names
colnames(pca.df.t6)[1:2] <- c("PC1", "PC2")

pdf("RAnalysis/Output/Time6.PCA.pdf", width=8, height=4)
par(mfrow=c(1,2))
# All Group Ordination
T6.PCA.A.plot <- ordiplot(T6.pca.out, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5), xlab=c("PC1 (42%)"), ylab=c("PC2 (27%)"))
points(T6.PCA.A.plot, "sites", col = groupcolors[pca.df.t6$group], cex=0.8, pch=PCA.symbols[pca.df.t6$Species])
ordihull(T6.PCA.A.plot, groups = pca.df.t6$group, draw = "polygon", alpha = 100, col = groupcolors, border = groupcolors)
par.new = T
plot(vect6, col = "black")
title(main = "A) Time 6, All Groups ")


#Transplant Ordination
T6.PCA.T.plot <- ordiplot(T6.pca.out, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5),xlab=c("PC1 (42%)"), ylab=c("PC2 (27%)"))
points(T6.PCA.T.plot, "sites", col = groupcolors[pca.df.t6$group], cex=0.8, pch=PCA.symbols[pca.df.t6$Species])
ordihull(T6.PCA.T.plot, groups = pca.df.t6$spe.tra, draw = "polygon", alpha = 100, col = PCA.colors, border = PCA.colors)
par.new = T
plot(vect6, col = "black")
title(main = "B) Time 6, Transplant Groups")
legend("bottomright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
                                 "PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
       pch = c(16,16,16,16,17,17,17,17),
       col = groupcolors, cex=0.4, box.lwd = 0,box.col = "white",bg = "white")
dev.off()


####### plasticity analysis #######
## T3
#plasticity (amount of change in multivariate space) is the distance between pairs of samples from the distance matrix of the response variables
T3.dist <- vegdist(t3scaled, method="euclidian")
#then need to extract distance for matched pairs 
#e.g., 81 vs 82 = parent colony 121 at the 2 different environments
T3.dist.df <- melt(as.matrix(T3.dist), varnames = c("identifier", "col"))
str(T3.dist.df)
T3.dist.df <- T3.dist.df[seq(2, nrow(T3.dist.df), 162), ]
T3.dist.df <- merge(T3.dist.df,t3, by="identifier")

#view summary stats
range(T3.dist.df$value)
boxplot(T3.dist.df$value)

#identify highest value
which(T3.dist.df$value ==max(T3.dist.df$value))

#plot genotype plasticity distance
par(mar=c(14,2,2 ,2))
boxplot(value ~group, data=T3.dist.df, las=2)
unique(T3.dist.df$group)
anova(aov(value ~group, data=T3.dist.df))

## T6
# one sample is missing, so need to deal with that in the coding
#plasticity (amount of change in multivariate space) is the distance between pairs of samples from the distance matrix of the response variables
T6.dist <- vegdist(t6scaled, method="euclidian")
#then need to extract distance for matched pairs 
#e.g., 81 vs 82 = parent colony 121 at the 2 different environments
T6.dist.df <- melt(as.matrix(T6.dist), varnames = c("identifier", "col"))
str(T6.dist.df)
T6.dist.df <- T6.dist.df[seq(2, nrow(T6.dist.df), 158), ]
T6.dist.df <- merge(T6.dist.df,t6, by="identifier")

#view summary stats
range(T6.dist.df$value)
boxplot(T6.dist.df$value)

#identify highest value
which(T6.dist.df$value ==max(T6.dist.df$value))

#plot genotype plasticity distance
par(mar=c(14,2,2 ,2))
boxplot(value ~group, data=T6.dist.df, las=2)
unique(T6.dist.df$group)
anova(aov(value ~group, data=T6.dist.df))

##### NMDS #####

# NMDS
## t3
nmds3 <- metaMDS(t3scaled, distance = "euclidian", autotransform = FALSE) #nmds
nmds3$stress
dev.off()
stressplot(nmds3)

## t6
nmds6 <- metaMDS(t6scaled, distance = "euclidian", autotransform = FALSE) #nmds
nmds6$stress
stressplot(nmds6)

# PERMANOVA
## t3
T3.mod <- adonis2(t3scaled ~ Transplant * Origin * Species, data = t3, method = "euclidian") # PERMANOVA
T3.mod
vect3 <- envfit(nmds3, t3[9:15], perm = 1000) #fit physiological vectors onto ordination
vect3df <- as.data.frame(vect3$vectors$arrows * sqrt(vect3$vectors$r))
vect3df$variable <- rownames(vect3df)

#set colors in this order MCAP Inner Lagoon MCAP Outer Lagoon PCOM Inner Lagoon PCOM Outer Lagoon
groupcolors <- c("cyan", "lightgreen", "green", "darkgreen", 
                 "pink", "purple", "orchid1", "darkorchid4") #set 8 group colors
MDS.colors <- c("lightgreen", "darkgreen", "orchid1", "darkorchid4") #set 4 Origin or Transplant colors
MDS.symbols <- c(16,17) #set two species shapes 16 = Mcap, 17 = Pcomp

nmds.df.t3 <- data.frame(nmds3$points[,1], y = nmds3$points[,2],
                         Parent = as.factor(fac.t3$Parent),
                         Origin = as.factor(fac.t3$Origin),
                         Transplant = as.factor(fac.t3$Transplant),
                         Species = as.factor(fac.t3$Species),
                         spe.tra = as.factor(fac.t3$spe.tra),
                         spe.ori = as.factor(fac.t3$spe.ori),
                         group = as.factor(fac.t3$group))

nmds.df.t3$origins <- factor(nmds.df.t3$spe.ori, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
nmds.df.t3$transplants <- factor(nmds.df.t3$spe.tra, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
nmds.df.t3$group <- factor(nmds.df.t3$group, levels = c("MCAP Inner Lagoon Inner Lagoon", "MCAP Inner Lagoon Outer Lagoon", "MCAP Outer Lagoon Outer Lagoon","MCAP Outer Lagoon Inner Lagoon", 
                                                        "PCOM Inner Lagoon Inner Lagoon", "PCOM Inner Lagoon Outer Lagoon", "PCOM Outer Lagoon Outer Lagoon", "PCOM Outer Lagoon Inner Lagoon"))
colnames(nmds.df.t3)[1:2] <- c("MDS1", "MDS2")

pdf("RAnalysis/Output/Time3.NMDS.pdf", width=8, height=8)
par(mfrow=c(2,2))
# All Group Ordination
#png("nmdst3vectors.All.png", height = 100, width = 150, units = "mm", res = 500)
T3.NMDS.plot <- ordiplot(nmds3, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T3.NMDS.plot, "sites", col = groupcolors[nmds.df.t3$group], cex=0.8, pch=MDS.symbols[nmds.df.t3$Species])
ordihull(T3.NMDS.plot, groups = nmds.df.t3$group, draw = "polygon", alpha = 100, col = groupcolors, border = groupcolors)
#ordiarrows(T3.NMDS.plot, groups = nmds.df.t3$Parent,  startmark = 1, label = nmds.df.t63Parent, length = .1)
par.new = T
plot(vect3, col = "black")
#legend("topright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
                                 #"PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
       #pch = c(16,16,16,16,17,17,17,17),
       #col = groupcolors, cex=0.6)
title(main = "A) Time 3, All Group Ordination")
#dev.off()

#Transplant Ordination
#png("nmdst3vectors.Transplant.png", height = 100, width = 150, units = "mm", res = 500)
T3.NMDS.plot <- ordiplot(nmds3, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T3.NMDS.plot, "sites", col = groupcolors[nmds.df.t3$group], cex=0.8, pch=MDS.symbols[nmds.df.t3$Species])
ordihull(T3.NMDS.plot, groups = nmds.df.t3$spe.tra, draw = "polygon", alpha = 100, col = MDS.colors, border = MDS.colors)
#ordiarrows(T3.NMDS.plot, groups = nmds.df.t3$Parent,  startmark = 1, label = nmds.df.t63Parent, length = .1)
par.new = T
plot(vect3, col = "black")
#legend("topright", legend = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"), 
       #pch = c(16,16,17,17),
       #col = MDS.colors, cex=0.6)
title(main = "B) Time 3, Transplant Ordination")
#dev.off()

# Group BetaDispersion Plasticity
#png("nmdst3vectors.Plasticity.png", height = 100, width = 150, units = "mm", res = 500)
T3.NMDS.plot <- ordiplot(nmds3, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T3.NMDS.plot, "sites", col = groupcolors[nmds.df.t3$group], cex=0.8, pch=MDS.symbols[nmds.df.t3$Species])
ordispider(nmds3, groups = nmds.df.t3$transplants, show.groups = nmds.df.t3$transplants, col = MDS.colors, label = F, spiders = c("centroid"))
#ordihull(T3.NMDS.plot, groups = nmds.df.t3$spe.tra, draw = "polygon", alpha = 100, col = MDS.colors, border = MDS.colors)
#ordiarrows(T3.NMDS.plot, groups = nmds.df.t3$Parent, col = "grey25", alpha = 50, length = .08)
#par.new = T
#plot(vect3, col = "black")
#legend("topright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
                              #"PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
       #pch = c(16,16,16,16,17,17,17,17),
       #col = groupcolors, cex=0.6)
title(main = "C) Time 3, Transplant Group Betadispersion")
#dev.off()


# Genotype Plasticity
#png("nmdst3vectors.Plasticity.png", height = 100, width = 150, units = "mm", res = 500)
T3.NMDS.plot <- ordiplot(nmds3, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T3.NMDS.plot, "sites", col = groupcolors[nmds.df.t3$group], cex=0.8, pch=MDS.symbols[nmds.df.t3$Species])
#ordihull(T3.NMDS.plot, groups = nmds.df.t3$spe.tra, draw = "polygon", alpha = 100, col = MDS.colors, border = MDS.colors)
ordiarrows(T3.NMDS.plot, groups = nmds.df.t3$Parent, col = "grey25", alpha = 50, length = .08)
#par.new = T
#plot(vect3, col = "black")
legend("bottomright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
                              "PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
       pch = c(16,16,16,16,17,17,17,17),
       col = groupcolors, cex=0.6, box.lwd = 0,box.col = "white",bg = "white")
title(main = "D) Time 3, Genotype plasticity")
dev.off()


bdisp <- betadisper(vegdist(t3scaled, method = "euclidian"), fac.t3$group)
anova(bdisp)


## t6
T6.mod <- adonis2(t6scaled ~ Transplant * Origin * Species, data = t6, method = "euclidian") # PERMANOVA
T6.mod
vect6 <- envfit(nmds6, t6[9:15], perm = 1000) #fit physiological vectors onto ordination
vect6df <- as.data.frame(vect6$vectors$arrows * sqrt(vect6$vectors$r))
vect6df$variable <- rownames(vect6df)

nmds.df.t6 <- data.frame(nmds6$points[,1], y = nmds6$points[,2],
                         Parent = as.factor(fac.t6$Parent),
                         Origin = as.factor(fac.t6$Origin),
                         Transplant = as.factor(fac.t6$Transplant),
                         Species = as.factor(fac.t6$Species),
                         spe.tra = as.factor(fac.t6$spe.tra),
                         spe.ori = as.factor(fac.t6$spe.ori),
                         group = as.factor(fac.t6$group))

nmds.df.t6$origins <- factor(nmds.df.t6$spe.ori, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
nmds.df.t6$transplants <- factor(nmds.df.t6$spe.tra, levels = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"))
nmds.df.t6$group <- factor(nmds.df.t6$group, levels = c("MCAP Inner Lagoon Inner Lagoon", "MCAP Inner Lagoon Outer Lagoon", "MCAP Outer Lagoon Outer Lagoon","MCAP Outer Lagoon Inner Lagoon", 
                                                        "PCOM Inner Lagoon Inner Lagoon", "PCOM Inner Lagoon Outer Lagoon", "PCOM Outer Lagoon Outer Lagoon", "PCOM Outer Lagoon Inner Lagoon"))
colnames(nmds.df.t6)[1:2] <- c("MDS1", "MDS2")

pdf("RAnalysis/Output/Time6.NMDS.pdf", width=8, height=8)
par(mfrow=c(2,2))
# All Group Ordination
#png("nmdst6vectors.All.png", height = 100, width = 150, units = "mm", res = 500)
T6.NMDS.plot <- ordiplot(nmds6, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T6.NMDS.plot, "sites", col = groupcolors[nmds.df.t6$group], cex=0.8, pch=MDS.symbols[nmds.df.t6$Species])
ordihull(T6.NMDS.plot, groups = nmds.df.t6$group, draw = "polygon", alpha = 100, col = groupcolors, border = groupcolors)
#ordiarrows(T6.NMDS.plot, groups = nmds.df.t6$Parent,  startmark = 1, label = nmds.df.t63Parent, length = .1)
par.new = T
plot(vect6, col = "black")
#legend("topright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
#"PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
#pch = c(16,16,16,16,17,17,17,17),
#col = groupcolors, cex=0.6)
title(main = "A) Time 6, All Group Ordination")
#dev.off()

#Transplant Ordination
#png("nmdst6vectors.Transplant.png", height = 100, width = 150, units = "mm", res = 500)
T6.NMDS.plot <- ordiplot(nmds6, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T6.NMDS.plot, "sites", col = groupcolors[nmds.df.t6$group], cex=0.8, pch=MDS.symbols[nmds.df.t6$Species])
ordihull(T6.NMDS.plot, groups = nmds.df.t6$spe.tra, draw = "polygon", alpha = 100, col = MDS.colors, border = MDS.colors)
#ordiarrows(T6.NMDS.plot, groups = nmds.df.t6$Parent,  startmark = 1, label = nmds.df.t63Parent, length = .1)
par.new = T
plot(vect6, col = "black")
#legend("topright", legend = c("MCAP Inner Lagoon", "MCAP Outer Lagoon", "PCOM Inner Lagoon", "PCOM Outer Lagoon"), 
#pch = c(16,16,17,17),
#col = MDS.colors, cex=0.6)
title(main = "B) Time 6, Transplant Ordination")
#dev.off()

# Group BetaDispersion Plasticity
#png("nmdst6vectors.Plasticity.png", height = 100, width = 150, units = "mm", res = 500)
T6.NMDS.plot <- ordiplot(nmds6, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T6.NMDS.plot, "sites", col = groupcolors[nmds.df.t6$group], cex=0.8, pch=MDS.symbols[nmds.df.t6$Species])
ordispider(nmds6, groups = nmds.df.t6$transplants, show.groups = nmds.df.t6$transplants, col = MDS.colors, label = F, spiders = c("centroid"))
#ordihull(T6.NMDS.plot, groups = nmds.df.t6$spe.tra, draw = "polygon", alpha = 100, col = MDS.colors, border = MDS.colors)
#ordiarrows(T6.NMDS.plot, groups = nmds.df.t6$Parent, col = "grey25", alpha = 50, length = .08)
#par.new = T
#plot(vect6, col = "black")
#legend("topright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
#"PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
#pch = c(16,16,16,16,17,17,17,17),
#col = groupcolors, cex=0.6)
title(main = "C) Time 6, Transplant Group Betadispersion")
#dev.off()


# Genotype Plasticity
#png("nmdst6vectors.Plasticity.png", height = 100, width = 150, units = "mm", res = 500)
T6.NMDS.plot <- ordiplot(nmds6, type = "n", display = "sites", ylim=c(-5,5), xlim=c(-5,5))
points(T6.NMDS.plot, "sites", col = groupcolors[nmds.df.t6$group], cex=0.8, pch=MDS.symbols[nmds.df.t6$Species])
#ordihull(T6.NMDS.plot, groups = nmds.df.t6$spe.tra, draw = "polygon", alpha = 100, col = MDS.colors, border = MDS.colors)
ordiarrows(T6.NMDS.plot, groups = nmds.df.t6$Parent, col = "grey25", alpha = 50, length = .08)
#par.new = T
#plot(vect6, col = "black")
legend("bottomright", legend = c("MCAP Inner:Inner", "MCAP Inner:Outer ", "MCAP Outer:Outer ","MCAP Outer:Inner ", 
                                 "PCOM Inner:Inner ", "PCOM Inner:Outer ", "PCOM Outer:Outer ", "PCOM Outer:Inner "), 
       pch = c(16,16,16,16,17,17,17,17),
       col = groupcolors, cex=0.6, box.lwd = 0,box.col = "white",bg = "white")
title(main = "D) Time 6, Genotype plasticity")
dev.off()


