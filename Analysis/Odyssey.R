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

rm(list = ls())

### Li-Cor Calibration ###
# Import files in a list
calibration.files <- list.files(path="data", pattern = "CSV$", full.names = T)

# Remove extensions from file names
file.names<-file_path_sans_ext(list.files(path="data", pattern = "CSV$", full.names = F))

# Formatting all data for loop
for(i in 1:length(calibration.files))
{
  data<-read.csv(calibration.files[i], sep=",", skip=8)
  df<-data[, c(-1,-5)]
  colnames(df)<-c("Date", "Time", "Raw.value")
  df$Date<-parse_date_time(df$Date, "dmy")
  df$timestamp<-as.POSIXct(paste(df$Date, df$Time), format="%Y-%m-%d %H:%M:%S")
  df<-df[!(df$timestamp < "2016-11-03 16:30:00"),] 
  df<-df[!(df$timestamp > "2016-11-04 13:00:00"),] 
  df<-df[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), df)) 
}

# List files and grab those that fill pattern 
files<-as.data.frame(mget(ls(pattern = "SN.*"))) 
data_index<-c(1,(seq(2,16,2))) 
all.data.PAR<-as.data.frame(c(files[, data_index])) 

# Name with logger SNs
colnames(all.data.PAR)<-c("timestamp", "SN2486", "SN6375", "SN6376", "SN6380", "SN6382", "SN6383", "SN6384", "SN7273")

# Add Licor file to the total logger dataframe 'all.data.PAR'
data<-read.csv("data/20161104_LicorCalibration.csv") 
df<-data[, c(4,5)] 
colnames(df)<-c("timestamp", "Licor.RAW")
df$timestamp<-mdy_hm(as.character(df$timestamp))
df$timestamp<-strptime(df$timestamp, format="%Y-%m-%d %H:%M:%S")

# Remove times where loggers were not logging with the LiCor
df<-df[!(df$timestamp < "2016-11-03 16:30:00"),] 
df<-df[!(df$timestamp > "2016-11-04 13:00:00"),]
df$Licor.RAW<-as.numeric(as.character(df$Licor.RAW)) 
Licor.data<-df 

# Trim to remove non-logged, dark periods according to Licor
all.data.PAR.sansnight<-all.data.PAR[!c(all.data.PAR$timestamp > "2016-11-03 16:30:00" & all.data.PAR$timestamp < "2016-11-04 08:15:00"), ]

# Integrate raw values over time interval (15min * 60 sec) and calibrate to LiCor
## Converts mol photons/interval to umol/photons/s
Licor.data$Licor.calibrated<-(Licor.data$Licor.RAW*(10^6)/(15*60))

# Final dataframe with Licor and Odyssey data
Ody.Licor.PAR<-as.data.frame(c(all.data.PAR.sansnight, Licor.data[3]))

# Generate standard curve with regression equation using "Ody.Licor.PAR" file
df<- Ody.Licor.PAR  
y<-Ody.Licor.PAR$Licor.calibrated
pdf(paste("plot", i,".pdf",sep=""))
par(mfrow=c(2,2))
for (i in 2:8) { 
  x<-Ody.Licor.PAR[,i] 
  m <- lm(y ~ x + 0, data = df) 
  plot(y ~ x, ylab="Cal. LiCor umol photons m-2 s-1", xlab="Odyssey reading", main=colnames(df)[i])
  abline(m)
  eq <-substitute(italic(y) == a + b~italic(x)*","~~italic(r)^2~"="~r2,
                  list(a = format(coef(m)[2], digits = 4),
                       b = format(coef(m)[1], digits = 4),
                       r2 = format(summary(m)$r.squared, digits = 3)))
  legend("topleft", legend=eq, bty="n")
}
dev.off()

###### Odyssey Data ######
# Import files in a list
R4F1_1 <- list.files(path="F1", pattern = "1006", full.names = T)
R4F1_2 <- list.files(path="F1", pattern = "1020", full.names = T)
R4F1_3 <- list.files(path="F1", pattern = "1202", full.names = T)
R4F1_4 <- list.files(path="F1", pattern = "026", full.names = T)
R4F1_5 <- list.files(path="F1", pattern = "027", full.names = T)
R4F1_6 <- list.files(path="F1", pattern = "20170315", full.names = T)
R4F1_7 <- list.files(path="F1", pattern = "March", full.names = T)
R4F2_1 <- list.files(path="F2", pattern = "1006", full.names = T)
R4F2_2 <- list.files(path="F2", pattern = "1020", full.names = T)
R4F2_3 <- list.files(path="F2", pattern = "1202", full.names = T)
R4F2_4 <- list.files(path="F2", pattern = "r4f24", full.names = T)
R4F2_5 <- list.files(path="F2", pattern = "frame", full.names = T)
R4F2_6 <- list.files(path="F2", pattern = "20170317", full.names = T)
R4F2_7 <- list.files(path="F2", pattern = "March", full.names = T)
R4F3_1 <- list.files(path="F3", pattern = "1006", full.names = T)
R4F3_2 <- list.files(path="F3", pattern = "1020", full.names = T)
R4F3_3 <- list.files(path="F3", pattern = "1202", full.names = T)
R13F4_1 <- list.files(path="F4", pattern = "1006", full.names = T)
R13F4_2 <- list.files(path="F4", pattern = "1020", full.names = T)
R13F4_3 <- list.files(path="F4", pattern = "1202", full.names = T)
R13F4_4 <- list.files(path="F4", pattern = "027", full.names = T)
R13F4_5 <- list.files(path="F4", pattern = "029", full.names = T)
R13F4_6 <- list.files(path="F4", pattern = "030", full.names = T)
R13F4_7 <- list.files(path="F4", pattern = "March", full.names = T)
R13F5_1 <- list.files(path="F5", pattern = "1006", full.names = T)
R13F5_2 <- list.files(path="F5", pattern = "1020", full.names = T)
R13F5_3 <- list.files(path="F5", pattern = "1202", full.names = T)
R13F5_4 <- list.files(path="F5", pattern = "r13f54", full.names = T)
R13F5_5 <- list.files(path="F5", pattern = "frame", full.names = T)
R13F6_1 <- list.files(path="F6", pattern = "1006", full.names = T)
R13F6_2 <- list.files(path="F6", pattern = "1020", full.names = T)
R13F6_3 <- list.files(path="F6", pattern = "1202", full.names = T)

# Remove extensions from file names
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "1006", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "1020", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "1202", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "026", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "027", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "20170315", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F1", pattern = "March", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "1006", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "1020", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "1202", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "r4f24", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "frame", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "20170317", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F2", pattern = "March", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F3", pattern = "1006", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F3", pattern = "1020", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F3", pattern = "1202", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "1006", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "1020", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "1202", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "027", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "029", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "030", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F4", pattern = "March", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F5", pattern = "1006", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F5", pattern = "1020", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F5", pattern = "1202", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F5", pattern = "r13f54", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F5", pattern = "frame", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F6", pattern = "1006", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F6", pattern = "1020", full.names = F))
file.names<-file_path_sans_ext(list.files(path="F6", pattern = "1202", full.names = F))

# Format all data - keep only timestamp and raw output
for(i in 1:length(R4F1_1))
{
  data<-read.csv(R4F1_1[i], sep=",", skip=8)
  F1_1<-data[, c(-1,-5)]
  colnames(F1_1)<-c("Date", "Time", "Raw.value")
  F1_1$Date<-parse_date_time(F1_1$Date, "dmy")
  F1_1$timestamp<-as.POSIXct(paste(F1_1$Date, F1_1$Time), format="%Y-%m-%d %H:%M:%S")
  F1_1<-F1_1[!(F1_1$timestamp < "2016-09-30 12:00:00"),] 
  F1_1<-F1_1[!(F1_1$timestamp > "2016-10-06 13:00:00"),] 
  F1_1<-F1_1[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_1)) 
}
names(F1_1) <- c("timestamp", "raw")

for(i in 1:length(R4F1_2))
{
  data<-read.csv(R4F1_2[i], sep=",", skip=8)
  F1_2<-data[, c(-1,-5)]
  colnames(F1_2)<-c("Date", "Time", "Raw.value")
  F1_2$Date<-parse_date_time(F1_2$Date, "dmy")
  F1_2$timestamp<-as.POSIXct(paste(F1_2$Date, F1_2$Time), format="%Y-%m-%d %H:%M:%S")
  F1_2<-F1_2[!(F1_2$timestamp < "2016-10-07 13:00:00"),] 
  F1_2<-F1_2[!(F1_2$timestamp > "2016-10-21 12:00:00"),] 
  F1_2<-F1_2[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_2)) 
}
names(F1_2) <- c("timestamp", "raw")

for(i in 1:length(R4F1_3))
{
  data<-read.csv(R4F1_3[i], sep=",", skip=8)
  F1_3<-data[, c(-1,-5)]
  colnames(F1_3)<-c("Date", "Time", "Raw.value")
  F1_3$Date<-parse_date_time(F1_3$Date, "dmy")
  F1_3$timestamp<-as.POSIXct(paste(F1_3$Date, F1_3$Time), format="%Y-%m-%d %H:%M:%S")
  F1_3<-F1_3[!(F1_3$timestamp < "2016-11-09 13:00:00"),] 
  F1_3<-F1_3[!(F1_3$timestamp > "2016-12-02 10:00:00"),] 
  F1_3<-F1_3[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_3)) 
}
names(F1_3) <- c("timestamp", "raw")

for(i in 1:length(R4F1_4))
{
  data<-read.csv(R4F1_4[i], sep=",", skip=8)
  F1_4<-data[, c(-1,-5)]
  colnames(F1_4)<-c("Date", "Time", "Raw.value")
  F1_4$Date<-parse_date_time(F1_4$Date, "dmy")
  F1_4$timestamp<-as.POSIXct(paste(F1_4$Date, F1_4$Time), format="%Y-%m-%d %H:%M:%S")
  F1_4<-F1_4[!(F1_4$timestamp < "2016-12-06 16:00:00"),] 
  F1_4<-F1_4[!(F1_4$timestamp > "2017-01-11 15:00:00"),] 
  F1_4<-F1_4[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_4)) 
}
names(F1_4) <- c("timestamp", "raw")

for(i in 1:length(R4F1_5))
{
  data<-read.csv(R4F1_5[i], sep=",", skip=8)
  F1_5<-data[, c(-1,-5)]
  colnames(F1_5)<-c("Date", "Time", "Raw.value")
  F1_5$Date<-parse_date_time(F1_5$Date, "dmy")
  F1_5$timestamp<-as.POSIXct(paste(F1_5$Date, F1_5$Time), format="%Y-%m-%d %H:%M:%S")
  F1_5<-F1_5[!(F1_5$timestamp < "2017-01-12 17:00:00"),] 
  F1_5<-F1_5[!(F1_5$timestamp > "2017-02-09 23:45:00"),] 
  F1_5<-F1_5[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_5)) 
}
names(F1_5) <- c("timestamp", "raw")

for(i in 1:length(R4F1_6))
{
  data<-read.csv(R4F1_6[i], sep=",", skip=8)
  F1_6<-data[, c(-1,-5)]
  colnames(F1_6)<-c("Date", "Time", "Raw.value")
  F1_6$Date<-parse_date_time(F1_6$Date, "dmy")
  F1_6$timestamp<-as.POSIXct(paste(F1_6$Date, F1_6$Time), format="%Y-%m-%d %H:%M:%S")
  F1_6<-F1_6[!(F1_6$timestamp < "2017-02-10 15:15:00"),] 
  F1_6<-F1_6[!(F1_6$timestamp > "2017-03-06 15:00:00"),] 
  F1_6<-F1_6[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_6)) 
}
names(F1_6) <- c("timestamp", "raw")

for(i in 1:length(R4F1_7))
{
  data<-read.csv(R4F1_7[i], sep=",", skip=8)
  F1_7<-data[, c(-1,-5)]
  colnames(F1_7)<-c("Date", "Time", "Raw.value")
  F1_7$Date<-parse_date_time(F1_7$Date, "dmy")
  F1_7$timestamp<-as.POSIXct(paste(F1_7$Date, F1_7$Time), format="%Y-%m-%d %H:%M:%S")
  F1_7<-F1_7[!(F1_7$timestamp < "2017-04-07 17:15:00"),] 
  F1_7<-F1_7[!(F1_7$timestamp > "2017-07-20 16:30:00"),] 
  F1_7<-F1_7[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F1_7)) 
}
names(F1_7) <- c("timestamp", "raw")

for(i in 1:length(R4F2_1))
{
  data<-read.csv(R4F2_1[i], sep=",", skip=9)
  F2_1<-data[, c(-1,-5)]
  colnames(F2_1)<-c("Date", "Time", "Raw.value")
  F2_1$Date<-parse_date_time(F2_1$Date, "dmy")
  F2_1$timestamp<-as.POSIXct(paste(F2_1$Date, F2_1$Time), format="%Y-%m-%d %H:%M:%S")
  F2_1<-F2_1[!(F2_1$timestamp < "2016-09-30 12:00:00"),] 
  F2_1<-F2_1[!(F2_1$timestamp > "2016-10-06 13:00:00"),] 
  F2_1<-F2_1[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_1)) 
}
names(F2_1) <- c("timestamp", "raw")

for(i in 1:length(R4F2_2))
{
  data<-read.csv(R4F2_2[i], sep=",", skip=8)
  F2_2<-data[, c(-1,-5)]
  colnames(F2_2)<-c("Date", "Time", "Raw.value")
  F2_2$Date<-parse_date_time(F2_2$Date, "dmy")
  F2_2$timestamp<-as.POSIXct(paste(F2_2$Date, F2_2$Time), format="%Y-%m-%d %H:%M:%S")
  F2_2<-F2_2[!(F2_2$timestamp < "2016-10-07 13:00:00"),] 
  F2_2<-F2_2[!(F2_2$timestamp > "2016-10-21 12:00:00"),] 
  F2_2<-F2_2[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_2)) 
}
names(F2_2) <- c("timestamp", "raw")

for(i in 1:length(R4F2_3))
{
  data<-read.csv(R4F2_3[i], sep=",", skip=8)
  F2_3<-data[, c(-1,-5)]
  colnames(F2_3)<-c("Date", "Time", "Raw.value")
  F2_3$Date<-parse_date_time(F2_3$Date, "dmy")
  F2_3$timestamp<-as.POSIXct(paste(F2_3$Date, F2_3$Time), format="%Y-%m-%d %H:%M:%S")
  F2_3<-F2_3[!(F2_3$timestamp < "2016-11-09 13:00:00"),] 
  F2_3<-F2_3[!(F2_3$timestamp > "2016-12-02 10:00:00"),] 
  F2_3<-F2_3[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_3)) 
}
names(F2_3) <- c("timestamp", "raw")

for(i in 1:length(R4F2_4))
{
  data<-read.csv(R4F2_4[i], sep=",", skip=8)
  F2_4<-data[, c(-1,-5)]
  colnames(F2_4)<-c("Date", "Time", "Raw.value")
  F2_4$Date<-parse_date_time(F2_4$Date, "dmy")
  F2_4$timestamp<-as.POSIXct(paste(F2_4$Date, F2_4$Time), format="%Y-%m-%d %H:%M:%S")
  F2_4<-F2_4[!(F2_4$timestamp < "2016-12-06 16:00:00"),] 
  F2_4<-F2_4[!(F2_4$timestamp > "2017-01-11 15:00:00"),] 
  F2_4<-F2_4[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_4)) 
}
names(F2_4) <- c("timestamp", "raw")

for(i in 1:length(R4F2_5))
{
  data<-read.csv(R4F2_5[i], sep=",", skip=8)
  F2_5<-data[, c(-1,-5)]
  colnames(F2_5)<-c("Date", "Time", "Raw.value")
  F2_5$Date<-parse_date_time(F2_5$Date, "dmy")
  F2_5$timestamp<-as.POSIXct(paste(F2_5$Date, F2_5$Time), format="%Y-%m-%d %H:%M:%S")
  F2_5<-F2_5[!(F2_5$timestamp < "2017-01-12 17:00:00"),] 
  F2_5<-F2_5[!(F2_5$timestamp > "2017-02-09 23:45:00"),] 
  F2_5<-F2_5[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_5)) 
}
names(F2_5) <- c("timestamp", "raw")

for(i in 1:length(R4F2_6))
{
  data<-read.csv(R4F2_6[i], sep=",", skip=8)
  F2_6<-data[, c(-1,-5)]
  colnames(F2_6)<-c("Date", "Time", "Raw.value")
  F2_6$Date<-parse_date_time(F2_6$Date, "dmy")
  F2_6$timestamp<-as.POSIXct(paste(F2_6$Date, F2_6$Time), format="%Y-%m-%d %H:%M:%S")
  F2_6<-F2_6[!(F2_6$timestamp < "2017-02-10 15:15:00"),] 
  F2_6<-F2_6[!(F2_6$timestamp > "2017-03-06 15:00:00"),] 
  F2_6<-F2_6[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_6)) 
}
names(F2_6) <- c("timestamp", "raw")

for(i in 1:length(R4F2_7))
{
  data<-read.csv(R4F2_7[i], sep=",", skip=8)
  F2_7<-data[, c(-1,-5)]
  colnames(F2_7)<-c("Date", "Time", "Raw.value")
  F2_7$Date<-parse_date_time(F2_7$Date, "dmy")
  F2_7$timestamp<-as.POSIXct(paste(F2_7$Date, F2_7$Time), format="%Y-%m-%d %H:%M:%S")
  F2_7<-F2_7[!(F2_7$timestamp < "2017-03-17 12:45:00"),] 
  F2_7<-F2_7[!(F2_7$timestamp > "2017-07-20 16:30:00"),] 
  F2_7<-F2_7[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F2_7)) 
}
names(F2_7) <- c("timestamp", "raw")

for(i in 1:length(R4F3_1))
{
  data<-read.csv(R4F3_1[i], sep=",", skip=9)
  F3_1<-data[, c(-1,-5)]
  colnames(F3_1)<-c("Date", "Time", "Raw.value")
  F3_1$Date<-parse_date_time(F3_1$Date, "dmy")
  F3_1$timestamp<-as.POSIXct(paste(F3_1$Date, F3_1$Time), format="%Y-%m-%d %H:%M:%S")
  F3_1<-F3_1[!(F3_1$timestamp < "2016-09-30 12:00:00"),] 
  F3_1<-F3_1[!(F3_1$timestamp > "2016-10-06 13:00:00"),] 
  F3_1<-F3_1[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F3_1)) 
}
names(F3_1) <- c("timestamp", "raw")

for(i in 1:length(R4F3_2))
{
  data<-read.csv(R4F3_2[i], sep=",", skip=8)
  F3_2<-data[, c(-1,-5)]
  colnames(F3_2)<-c("Date", "Time", "Raw.value")
  F3_2$Date<-parse_date_time(F3_2$Date, "dmy")
  F3_2$timestamp<-as.POSIXct(paste(F3_2$Date, F3_2$Time), format="%Y-%m-%d %H:%M:%S")
  F3_2<-F3_2[!(F3_2$timestamp < "2016-10-07 13:00:00"),] 
  F3_2<-F3_2[!(F3_2$timestamp > "2016-10-21 12:00:00"),] 
  F3_2<-F3_2[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F3_2)) 
}
names(F3_2) <- c("timestamp", "raw")

for(i in 1:length(R4F3_3))
{
  data<-read.csv(R4F3_3[i], sep=",", skip=8)
  F3_3<-data[, c(-1,-5)]
  colnames(F3_3)<-c("Date", "Time", "Raw.value")
  F3_3$Date<-parse_date_time(F3_3$Date, "dmy")
  F3_3$timestamp<-as.POSIXct(paste(F3_3$Date, F3_3$Time), format="%Y-%m-%d %H:%M:%S")
  F3_3<-F3_3[!(F3_3$timestamp < "2016-11-09 13:00:00"),] 
  F3_3<-F3_3[!(F3_3$timestamp > "2016-12-02 10:00:00"),] 
  F3_3<-F3_3[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F3_3)) 
}
names(F3_3) <- c("timestamp", "raw")

for(i in 1:length(R13F4_1))
{
  data<-read.csv(R13F4_1[i], sep=",", skip=8)
  F4_1<-data[, c(-1,-5)]
  colnames(F4_1)<-c("Date", "Time", "Raw.value")
  F4_1$Date<-parse_date_time(F4_1$Date, "dmy")
  F4_1$timestamp<-as.POSIXct(paste(F4_1$Date, F4_1$Time), format="%Y-%m-%d %H:%M:%S")
  F4_1<-F4_1[!(F4_1$timestamp < "2016-09-30 13:00:00"),] 
  F4_1<-F4_1[!(F4_1$timestamp > "2016-10-06 14:30:00"),] 
  F4_1<-F4_1[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_1)) 
}
names(F4_1) <- c("timestamp", "raw")

for(i in 1:length(R13F4_2))
{
  data<-read.csv(R13F4_2[i], sep=",", skip=8)
  F4_2<-data[, c(-1,-5)]
  colnames(F4_2)<-c("Date", "Time", "Raw.value")
  F4_2$Date<-parse_date_time(F4_2$Date, "dmy")
  F4_2$timestamp<-as.POSIXct(paste(F4_2$Date, F4_2$Time), format="%Y-%m-%d %H:%M:%S")
  F4_2<-F4_2[!(F4_2$timestamp < "2016-10-07 14:00:00"),] 
  F4_2<-F4_2[!(F4_2$timestamp > "2016-10-21 12:00:00"),] 
  F4_2<-F4_2[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_2)) 
}
names(F4_2) <- c("timestamp", "raw")

for(i in 1:length(R13F4_3))
{
  data<-read.csv(R13F4_3[i], sep=",", skip=8)
  F4_3<-data[, c(-1,-5)]
  colnames(F4_3)<-c("Date", "Time", "Raw.value")
  F4_3$Date<-parse_date_time(F4_3$Date, "dmy")
  F4_3$timestamp<-as.POSIXct(paste(F4_3$Date, F4_3$Time), format="%Y-%m-%d %H:%M:%S")
  F4_3<-F4_3[!(F4_3$timestamp < "2016-11-09 14:00:00"),] 
  F4_3<-F4_3[!(F4_3$timestamp > "2016-12-02 11:00:00"),] 
  F4_3<-F4_3[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_3)) 
}
names(F4_3) <- c("timestamp", "raw")

for(i in 1:length(R13F4_4))
{
  data<-read.csv(R13F4_4[i], sep=",", skip=8)
  F4_4<-data[, c(-1,-5)]
  colnames(F4_4)<-c("Date", "Time", "Raw.value")
  F4_4$Date<-parse_date_time(F4_4$Date, "dmy")
  F4_4$timestamp<-as.POSIXct(paste(F4_4$Date, F4_4$Time), format="%Y-%m-%d %H:%M:%S")
  F4_4<-F4_4[!(F4_4$timestamp < "2016-12-05 16:00:00"),] 
  F4_4<-F4_4[!(F4_4$timestamp > "2017-01-11 15:00:00"),] 
  F4_4<-F4_4[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_4)) 
}
names(F4_4) <- c("timestamp", "raw")

for(i in 1:length(R13F4_5))
{
  data<-read.csv(R13F4_5[i], sep=",", skip=8)
  F4_5<-data[, c(-1,-5)]
  colnames(F4_5)<-c("Date", "Time", "Raw.value")
  F4_5$Date<-parse_date_time(F4_5$Date, "dmy")
  F4_5$timestamp<-as.POSIXct(paste(F4_5$Date, F4_5$Time), format="%Y-%m-%d %H:%M:%S")
  F4_5<-F4_5[!(F4_5$timestamp < "2017-01-12 17:00:00"),] 
  F4_5<-F4_5[!(F4_5$timestamp > "2017-02-09 23:45:00"),] 
  F4_5<-F4_5[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_5)) 
}
names(F4_5) <- c("timestamp", "raw")

for(i in 1:length(R13F4_6))
{
  data<-read.csv(R13F4_6[i], sep=",", skip=8)
  F4_6<-data[, c(-1,-5)]
  colnames(F4_6)<-c("Date", "Time", "Raw.value")
  F4_6$Date<-parse_date_time(F4_6$Date, "dmy")
  F4_6$timestamp<-as.POSIXct(paste(F4_6$Date, F4_6$Time), format="%Y-%m-%d %H:%M:%S")
  F4_6<-F4_6[!(F4_6$timestamp < "2017-02-10 15:15:00"),] 
  F4_6<-F4_6[!(F4_6$timestamp > "2017-03-06 15:00:00"),] 
  F4_6<-F4_6[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_6)) 
}
names(F4_6) <- c("timestamp", "raw")

for(i in 1:length(R13F4_7))
{
  data<-read.csv(R13F4_7[i], sep=",", skip=8)
  F4_7<-data[, c(-1,-5:-10)]
  colnames(F4_7)<-c("Date", "Time", "Raw.value")
  F4_7$Date<-parse_date_time(F4_7$Date, "dmy")
  F4_7$timestamp<-as.POSIXct(paste(F4_7$Date, F4_7$Time), format="%Y-%m-%d %H:%M:%S")
  F4_7<-F4_7[!(F4_7$timestamp < "2017-03-17 12:45:00"),] 
  F4_7<-F4_7[!(F4_7$timestamp > "2017-07-20 16:45:00"),] 
  F4_7<-F4_7[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F4_7)) 
}
names(F4_7) <- c("timestamp", "raw")

for(i in 1:length(R13F5_1))
{
  data<-read.csv(R13F5_1[i], sep=",", skip=8)
  F5_1<-data[, c(-1,-5)]
  colnames(F5_1)<-c("Date", "Time", "Raw.value")
  F5_1$Date<-parse_date_time(F5_1$Date, "dmy")
  F5_1$timestamp<-as.POSIXct(paste(F5_1$Date, F5_1$Time), format="%Y-%m-%d %H:%M:%S")
  F5_1<-F5_1[!(F5_1$timestamp < "2016-09-30 13:00:00"),] 
  F5_1<-F5_1[!(F5_1$timestamp > "2016-10-06 14:30:00"),] 
  F5_1<-F5_1[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F5_1)) 
}
names(F5_1) <- c("timestamp", "raw")

for(i in 1:length(R13F5_2))
{
  data<-read.csv(R13F5_2[i], sep=",", skip=8)
  F5_2<-data[, c(-1,-5)]
  colnames(F5_2)<-c("Date", "Time", "Raw.value")
  F5_2$Date<-parse_date_time(F5_2$Date, "dmy")
  F5_2$timestamp<-as.POSIXct(paste(F5_2$Date, F5_2$Time), format="%Y-%m-%d %H:%M:%S")
  F5_2<-F5_2[!(F5_2$timestamp < "2016-10-07 14:00:00"),] 
  F5_2<-F5_2[!(F5_2$timestamp > "2016-10-21 12:00:00"),] 
  F5_2<-F5_2[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F5_2)) 
}
names(F5_2) <- c("timestamp", "raw")

for(i in 1:length(R13F5_3))
{
  data<-read.csv(R13F5_3[i], sep=",", skip=8)
  F5_3<-data[, c(-1,-5)]
  colnames(F5_3)<-c("Date", "Time", "Raw.value")
  F5_3$Date<-parse_date_time(F5_3$Date, "dmy")
  F5_3$timestamp<-as.POSIXct(paste(F5_3$Date, F5_3$Time), format="%Y-%m-%d %H:%M:%S")
  F5_3<-F5_3[!(F5_3$timestamp < "2016-11-09 14:00:00"),] 
  F5_3<-F5_3[!(F5_3$timestamp > "2016-12-02 11:00:00"),] 
  F5_3<-F5_3[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F5_3)) 
}
names(F5_3) <- c("timestamp", "raw")

for(i in 1:length(R13F5_4))
{
  data<-read.csv(R13F5_4[i], sep=",", skip=8)
  F5_4<-data[, c(-1,-5)]
  colnames(F5_4)<-c("Date", "Time", "Raw.value")
  F5_4$Date<-parse_date_time(F5_4$Date, "dmy")
  F5_4$timestamp<-as.POSIXct(paste(F5_4$Date, F5_4$Time), format="%Y-%m-%d %H:%M:%S")
  F5_4<-F5_4[!(F5_4$timestamp < "2016-12-05 17:00:00"),] 
  F5_4<-F5_4[!(F5_4$timestamp > "2017-01-11 15:00:00"),] 
  F5_4<-F5_4[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F5_4)) 
}
names(F5_4) <- c("timestamp", "raw")

for(i in 1:length(R13F5_5))
{
  data<-read.csv(R13F5_5[i], sep=",", skip=8)
  F5_5<-data[, c(-1,-5)]
  colnames(F5_5)<-c("Date", "Time", "Raw.value")
  F5_5$Date<-parse_date_time(F5_5$Date, "dmy")
  F5_5$timestamp<-as.POSIXct(paste(F5_5$Date, F5_5$Time), format="%Y-%m-%d %H:%M:%S")
  F5_5<-F5_5[!(F5_5$timestamp < "2017-01-12 17:00:00"),] 
  F5_5<-F5_5[!(F5_5$timestamp > "2017-02-09 23:45:00"),] 
  F5_5<-F5_5[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F5_5)) 
}
names(F5_5) <- c("timestamp", "raw")

for(i in 1:length(R13F6_1))
{
  data<-read.csv(R13F6_1[i], sep=",", skip=9)
  F6_1<-data[, c(-1,-5)]
  colnames(F6_1)<-c("Date", "Time", "Raw.value")
  F6_1$Date<-parse_date_time(F6_1$Date, "dmy")
  F6_1$timestamp<-as.POSIXct(paste(F6_1$Date, F6_1$Time), format="%Y-%m-%d %H:%M:%S")
  F6_1<-F6_1[!(F6_1$timestamp < "2016-09-30 12:00:00"),] 
  F6_1<-F6_1[!(F6_1$timestamp > "2016-10-06 13:00:00"),] 
  F6_1<-F6_1[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F6_1)) 
}
names(F6_1) <- c("timestamp", "raw")

for(i in 1:length(R13F6_2))
{
  data<-read.csv(R13F6_2[i], sep=",", skip=8)
  F6_2<-data[, c(-1,-5)]
  colnames(F6_2)<-c("Date", "Time", "Raw.value")
  F6_2$Date<-parse_date_time(F6_2$Date, "dmy")
  F6_2$timestamp<-as.POSIXct(paste(F6_2$Date, F6_2$Time), format="%Y-%m-%d %H:%M:%S")
  F6_2<-F6_2[!(F6_2$timestamp < "2016-10-07 13:00:00"),] 
  F6_2<-F6_2[!(F6_2$timestamp > "2016-10-21 12:00:00"),] 
  F6_2<-F6_2[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F6_2)) 
}
names(F6_2) <- c("timestamp", "raw")

for(i in 1:length(R13F6_3))
{
  data<-read.csv(R13F6_3[i], sep=",", skip=8)
  F6_3<-data[, c(-1,-5)]
  colnames(F6_3)<-c("Date", "Time", "Raw.value")
  F6_3$Date<-parse_date_time(F6_3$Date, "dmy")
  F6_3$timestamp<-as.POSIXct(paste(F6_3$Date, F6_3$Time), format="%Y-%m-%d %H:%M:%S")
  F6_3<-F6_3[!(F6_3$timestamp < "2016-11-09 13:00:00"),] 
  F6_3<-F6_3[!(F6_3$timestamp > "2016-12-02 10:00:00"),] 
  F6_3<-F6_3[ , c(4,3)] 
  make.names(assign(paste("SN",file.names[i], sep=""), F6_3)) 
}
names(F6_3) <- c("timestamp", "raw")

# Bind dataframes by reef
F1 <- rbind(F1_1, F1_2, F1_3, F1_4, F1_5, F1_6, F1_7)
F1$cal <- F1$raw*0.1317
F1$Rack <- "F1"
F1$Date <- as.Date(F1$timestamp) 
F1$Time <- format(as.POSIXct(F1$timestamp) ,format = "%H:%M:%S") 
F1NoNight <- with(F1, F1[hour(timestamp) > 4,])
F1NoNight <- with(F1NoNight, F1NoNight[hour(timestamp) < 20,])
F2 <- rbind(F2_1, F2_2, F2_3, F2_4, F2_5, F2_6, F2_7)
F2$cal <- F2$raw*0.119
F2$Rack <- "F2"
F2$Date <- as.Date(F2$timestamp) 
F2$Time <- format(as.POSIXct(F2$timestamp) ,format = "%H:%M:%S") 
F2NoNight <- with(F2, F2[hour(timestamp) > 4,])
F2NoNight <- with(F2NoNight, F2NoNight[hour(timestamp) < 20,])
F3 <- rbind(F3_1, F3_2, F3_3)
F3$cal <- F3$raw*0.1021
F3$Rack <- "F3"
F3$Date <- as.Date(F3$timestamp) 
F3$Time <- format(as.POSIXct(F3$timestamp) ,format = "%H:%M:%S") 
F3NoNight <- with(F3, F3[hour(timestamp) > 4,])
F3NoNight <- with(F3NoNight, F3NoNight[hour(timestamp) < 20,])
F4 <- rbind(F4_1, F4_2, F4_3, F4_4, F4_5, F4_6, F4_7)
F4$cal <- F4$raw*0.09455
F4$Rack <- "F4"
F4$Date <- as.Date(F4$timestamp) 
F4$Time <- format(as.POSIXct(F4$timestamp) ,format = "%H:%M:%S") 
F4NoNight <- with(F4, F4[hour(timestamp) > 4,])
F4NoNight <- with(F4NoNight, F4NoNight[hour(timestamp) < 20,])
F5 <- rbind(F5_1, F5_2, F5_3, F5_5)
F5$cal <- F5$raw*0.1193
F5$Rack <- "F5"
F5$Date <- as.Date(F5$timestamp) 
F5$Time <- format(as.POSIXct(F5$timestamp) ,format = "%H:%M:%S") 
F5NoNight <- with(F5, F5[hour(timestamp) > 4,])
F5NoNight <- with(F5NoNight, F5NoNight[hour(timestamp) < 20,])
F6 <- rbind(F6_1, F6_2, F6_3)
F6$cal <- F6$raw*0.09052
F6$Rack <- "F6"
F6$Date <- as.Date(F6$timestamp) 
F6$Time <- format(as.POSIXct(F6$timestamp) ,format = "%H:%M:%S") 
F6NoNight <- with(F6, F6[hour(timestamp) > 4,])
F6NoNight <- with(F6NoNight, F6NoNight[hour(timestamp) < 20,])

# Merge frames by site and get timestamp mean
sem <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))
ci <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))*qt(0.975, length(x)-1)
range <- function(x) max(x, na.rm = T)-min(x, na.rm = T)

PAR <- rbind(F1, F2, F3, F4, F5, F6)
PAR$Site <- ifelse(PAR$Rack == "F1" | PAR$Rack == "F2" | PAR$Rack == "F3", "IB", "OB")
MeanPAR <- cbind(aggregate(PAR$cal, by = list(PAR$Site, PAR$Time), FUN = mean), aggregate(PAR$cal, by = list(PAR$Site, PAR$Time), FUN = sem), aggregate(PAR$cal, by = list(PAR$Site, PAR$Time), FUN = ci))
MeanPAR <- MeanPAR[, c(1:3,6,9)]
names(MeanPAR) <- c("Site", "Time", "PAR", "PARSEM", "PARCI")
MeanPAR <- subset(MeanPAR, Time >= "05:00:00" & Time <= "20:00:00")
write.csv(PAR, "PAR.csv")

# Subset by timepoint
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

