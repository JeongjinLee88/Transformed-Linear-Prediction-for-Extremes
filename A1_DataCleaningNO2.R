#setwd("")
#setwd("/home/leej40/Documents/extlinear/code")

####  Bind multiple csv files
library(readr)
library(dplyr)
files <- list.files(path = "C:/NO2", pattern = "*.csv", full.names = T)
NO2 <- sapply(files, read_csv, simplify=FALSE) %>%
  bind_rows(.id = "id")

NO2=NO2[,c(2,4,6)]  # Date, Site ID, Max 1hr NO2
NO2=as.data.frame(NO2)
NO2$Date=as.Date(NO2$Date,"%m/%d/%Y")
str(NO2)

####  Extract each monitor separately
Alx=NO2[NO2[,2]==515100009,c(1,3)]
Mc=NO2[NO2[,2]==110010043,c(1,3)]
Rt=NO2[NO2[,2]==110010041,c(1,3)]
Ts=NO2[NO2[,2]==110010025,c(1,3)]
Trc=NO2[NO2[,2]==110010050,c(1,3)]
Arl=NO2[NO2[,2]==510130020,c(1,3)]
City_Alx=NO2[NO2[,2]==515100021,c(1,3)]

colnames(x = Alx)=c("Date","Max1hr_NO2_Alx")
colnames(x = Mc)=c("Date","Max1hr_NO2_Mc")
colnames(x = Rt)=c("Date","Max1hr_NO2_Rt")
colnames(x = Ts)=c("Date","Max1hr_NO2_Ts")
colnames(x = Trc)=c("Date","Max1hr_NO2_Ts")
colnames(x = Arl)=c("Date","Max1hr_NO2_Arl")
colnames(x = City_Alx)=c("Date","Max1hr_NO2_Alx")

####  Check the shape parameter
#library(evd)
library(ismev)
dev.new()
mrl.plot(Alx[,2])
q=quantile(Alx[,2],probs = c(0.9,0.95,0.99))
abline(v=q,col=c(1,2,3))
Alx_fit=gpd.fit(Alx[,2],threshold = 60)
#[1] 8.86191982 0.01074944
gpd.diag(z = Alx_fit)
gpd.profxi(Alx_fit,-0.1,0.5)

mrl.plot(City_Alx[,2])
q=quantile(City_Alx[,2],probs = c(0.9,0.95,0.99))
abline(v=q,col=c(1,2,3))
CityAlx_fit=gpd.fit(City_Alx[,2],threshold = 49.6)
#[1] 12.18869479 -0.09412186
gpd.diag(z = CityAlx_fit)
gpd.profxi(CityAlx_fit,-0.2,0.4)

mrl.plot(Mc[,2])
quantile(Mc[,2],probs = c(0.9,0.95,0.99))
abline(v=56)
Mc_fit=gpd.fit(Mc[,2],threshold = 56)
#[1] 7.9238372 0.1365265
gpd.diag(z = Mc_fit)
gpd.profxi(Mc_fit,-0.1,0.2)

mrl.plot(Rt[,2])
quantile(Rt[,2],probs = c(0.9,0.95,0.99))
abline(v=59)
Rt_fit=gpd.fit(Rt[,2],threshold = 59)
#[1] 9.81123356 0.07020839
gpd.diag(z = Rt_fit)
gpd.profxi(Rt_fit,-0.1,0.2)

####  Total measurements for each station
View(Alx[,1]) #1995-2012
View(City_Alx[,1]) #2012-2016
View(Mc[,1]) #1995-2020
View(Rt[,1]) #1995-2014,Mar / 2016, June - 2020,April
View(Ts[,1]) #1995-2010
View(Trc[,1]) #2013-2020
View(Arl[,1]) #1995-2020

c(dim(Alx)[1],dim(Mc)[1],dim(Rt)[1],dim(Ts)[1],dim(Trc)[1],dim(Arl)[1],dim(City_Alx)[1])
#[1] 6062 9033 8291 5725 2513 8838 1120

####  Combine Tacoma school with Tacoma rec center
TS_combined=rbind(x = Ts,y = Trc)
TS_combined$color=as.character(cut(TS_combined[,1],breaks = c(as.Date("1995-02-01"),as.Date("2010-12-31"),as.Date("2020-04-30")),labels = c("black","blue"),right = FALSE))
####  Merge McMill with Tacoma school+Tacoma rec center
Mc_Ts=merge(x = Mc,y = TS_combined,by = "Date")
####  Plot McMill vs Tacoma school+Tacoma rec center
dev.new()
par(mar=c(5.1,5.1,2,2))
plot(Mc_Ts[,2],Mc_Ts[,3],col=Mc_Ts$color,main="",xlab="McMillan",ylab="Tacoma school")
####  Combine Alexandria with City of Alexandria
Alx_combined=rbind(Alx,City_Alx)  #6062:
Alx_combined$color=as.character(cut(Alx_combined[,1],breaks = c(as.Date("1995-01-01"),as.Date("2012-08-20"),as.Date("2016-04-30")),labels = c("black","blue"),right = FALSE))
####  Merge McMill with Alx
Mc_Alx=merge(x = Mc,y = Alx_combined,by = "Date")
####  Plot McMill vs Tacoma school+Tacoma rec center
dev.new()
par(mar=c(5.1,5.1,2,2))
plot(Mc_Alx[,2],Mc_Alx[,3],col=Mc_Alx$color,main="",xlab="McMillan",ylab="Alexandria", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Old_Alx=Mc_Alx[Mc_Alx[,4]=="black",]
New_Alx=Mc_Alx[Mc_Alx[,4]=="blue",]
plot(Old_Alx[,2],Old_Alx[,3],ylim=c(0,140),xlim=c(0,255),col="black",main="",xlab="McMillan",ylab="Alexandria", cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=0.8,pch=1)
par(new=TRUE)
plot(New_Alx[,2],New_Alx[,3],ylim=c(0,140),xlim=c(0,255),col="blue",main="",xlab="McMillan",ylab="Alexandria", cex.main=1.5, cex.lab=1.5, cex.axis=1.5,pch=2)

Alx_2012=subset(Alx,Alx$Date >= as.Date("2012-01-01") & Alx$Date <= as.Date("2012-12-31"))
City_Alx2012=subset(City_Alx,City_Alx$Date >= as.Date("2012-01-01") & City_Alx$Date <= as.Date("2012-12-31"))
Alx_08=subset(Alx,Alx$Date >= as.Date("2012-08-01") & Alx$Date <= as.Date("2012-08-31"))
City_Alx08=subset(City_Alx,City_Alx$Date >= as.Date("2012-08-24") & City_Alx$Date <= as.Date("2012-08-31"))

dev.new()
par(mfrow=c(1,2))
plot(Alx_2012,type="l",ylim=c(0,80))
plot(City_Alx2012,type="l",ylim=c(0,80))

summary(Alx_2012[,2])
summary(City_Alx2012[,2])
summary(Alx_08[,2])
summary(City_Alx08[,2])
##  Dates are unique. There is no duplicate obs

####  Merge multiple data frames
NO2_fin=Reduce(function(x,y) merge(x = x, y = y, by = "Date",sort = T), 
               list(Alx_combined[,1:2], Mc, Rt, TS_combined[,1:2], Arl))

c(dim(Alx_combined)[1],dim(Mc)[1],dim(Rt)[1],dim(TS_combined)[1],dim(Arl)[1])
dim(NO2_fin)[1] #Total 5163

####  Standardize NO2 at Alexandria
source("A2_Mov_Avg.R")
####  Plot of standardized NO2 at Alx
dev.new()
par(mar=c(5.1,5.1,2,2))
plot(Alx_combined[,2],type="l",xlab="Date",ylab="NO2 measurements at Alexandria",ylim=c(0,150),xaxt = "n",col=grey(0.5), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Alx_MA=Mov_Avg(Alx_combined[,2], 901, TRUE)
lines(Alx_MA$MA_mean, col="blue", lwd=3,lty=2)
lines(Alx_MA$MA_sd, col="red", lwd=2)

####  Standardize NO2 at Alexandria via MA trend method
Std_Alx=(Alx_combined[,2]-Alx_MA$MA_mean)/Alx_MA$MA_sd
dev.new()
par(mar=c(5.1,5.1,2,2))
plot(Std_Alx,type="l",xlab="Date",ylab="Standardized NO2 measurements at Alexandria",xaxt="n",col=grey(0.5), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Alx_MA=Mov_Avg(Std_Alx, 901, TRUE)
lines(Alx_MA$MA_mean, col="blue", lwd=3,lty=2)
lines(Alx_MA$MA_sd, col="red", lwd=2)

save(Std_Alx,file="Std_Alx.Rdata")

####  Detect any trends
#library(forecast)
#library(tstools)
#library(devtools)

####  Time series plots
source("A2_Mov_Avg.R")
dev.new()
par(mar=c(5.1,5.1,2,2))
par(mfrow=c(1,2))
plot(NO2_fin[,2],type="l",xlab="Date",ylab="NO2 measurements at Alexandria",ylim=c(0,150),xaxt = "n",col=grey(0.5), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Alx_MA=Mov_Avg(NO2_fin[,2], 901, TRUE)
lines(Alx_MA$MA_mean, col="blue", lwd=3,lty=2)
lines(Alx_MA$MA_sd, col="red", lwd=2)

####  Check ACF and PACF
dev.new()
par(mfrow=c(1,2))
acf(NO2_fin[,2])
pacf(NO2_fin[,2])

####  Standardize NO2 at Alexandria
Std_Alx=(NO2_fin[,2]-Alx_MA$MA_mean)/Alx_MA$MA_sd
plot(Std_Alx,type="l",xlab="Date",ylab="Standardized NO2 measurements at Alexandria",xaxt="n",col=grey(0.5), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Alx_MA=Mov_Avg(Std_Alx, 901, TRUE)
lines(Alx_MA$MA_mean, col="blue", lwd=3,lty=2)
lines(Alx_MA$MA_sd, col="red", lwd=2)

####  Check ACF and PACF
dev.new()
par(mfrow=c(1,2))
acf(Std_Alx)
pacf(Std_Alx)

####  McMillan
dev.new()
plot(NO2_fin[,3],type="l",xlab="Year",ylab="NO2 measurements",ylim=c(0,150),xaxt="n")
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010))
Mc_MA=Mov_Avg(NO2_fin[,3], 901, TRUE)
lines(Mc_MA$MA_mean, col="green", lwd=2)
lines(Mc_MA$MA_sd, col="red", lwd=2)

dev.new()
par(mfrow=c(1,2))
acf(NO2_fin[,3])
pacf(NO2_fin[,3])

####  Standardize NO2 at McMillan
Std_Mc=(NO2_fin[,3]-Mc_MA$MA_mean)/Mc_MA$MA_sd
dev.new()
plot(Std_Mc,type="l",xlab="Year",ylab="NO2 measurements",xaxt="n")
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010))
Mc_MA=Mov_Avg(Std_Mc, 901, TRUE)
lines(Mc_MA$MA_mean, col="green", lwd=2)
lines(Mc_MA$MA_sd, col="red", lwd=2)

####  Check ACF and PACF
dev.new()
par(mfrow=c(1,2))
acf(Std_Mc)
pacf(Std_Mc)

####  River Terrace
dev.new()
plot(NO2_fin[,4],type="l",xlab="Year",ylab="NO2 measurements",ylim=c(0,150),xaxt="n")
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010))
Rt_MA=Mov_Avg(NO2_fin[,4], 901, TRUE)
lines(Rt_MA$MA_mean, col="green", lwd=2)
lines(Rt_MA$MA_sd, col="red", lwd=2)

dev.new()
par(mfrow=c(1,2))
acf(NO2_fin[,4])
pacf(NO2_fin[,4])

####  Standardize NO2 at River terrace
Std_Rt=(NO2_fin[,4]-Rt_MA$MA_mean)/Rt_MA$MA_sd
dev.new()
plot(Std_Rt,type="l",xlab="Year",ylab="NO2 measurements",xaxt="n")
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010))
Rt_MA=Mov_Avg(Std_Rt, 901, TRUE)
lines(Rt_MA$MA_mean, col="green", lwd=2)
lines(Rt_MA$MA_sd, col="red", lwd=2)

####  Check ACF and PACF
dev.new()
par(mfrow=c(1,2))
acf(Std_Rt)
pacf(Std_Rt)

####  Tacoma School
dev.new()
plot(NO2_fin[,5],type="l",xlab="Year",ylab="NO2 measurements",ylim=c(0,150),xaxt="n")
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010))
Ts_MA=Mov_Avg(NO2_fin[,5], 901, TRUE)
lines(Ts_MA$MA_mean, col="green", lwd=2)
lines(Ts_MA$MA_sd, col="red", lwd=2)

dev.new()
par(mfrow=c(1,2))
acf(NO2_fin[,5])
pacf(NO2_fin[,5])

####  Standardize NO2 at Tacoma school
Std_Ts=(NO2_fin[,5]-Ts_MA$MA_mean)/Ts_MA$MA_sd
dev.new()
plot(Std_Ts,type="l",xlab="Year",ylab="NO2 measurements",xaxt="n")
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010))
Ts_MA=Mov_Avg(Std_Ts, 901, TRUE)
lines(Ts_MA$MA_mean, col="green", lwd=2)
lines(Ts_MA$MA_sd, col="red", lwd=2)

####  Check ACF and PACF
dev.new()
par(mfrow=c(1,2))
acf(Std_Ts)
pacf(Std_Ts)

####  Arlignton
dev.new()
par(mar=c(5.1,5.1,2,2))
par(mfrow=c(1,2))
plot(NO2_fin[,6],type="l",col="grey40",xlab="Year",ylab="NO2 measurements at Arlington",ylim=c(0,150),xaxt="n", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Arl_MA=Mov_Avg(NO2_fin[,6], 901, TRUE)
lines(Arl_MA$MA_mean, col="blue", lwd=3,lty=2)
lines(Arl_MA$MA_sd, col="red", lwd=2)

dev.new()
par(mfrow=c(1,2))
acf(NO2_fin[,6])
pacf(NO2_fin[,6])

####  Standardize NO2 at Arlignton
Std_Arl=(NO2_fin[,6]-Arl_MA$MA_mean)/Arl_MA$MA_sd
plot(Std_Arl,type="l",col="grey40",xlab="Year",ylab="Standardized NO2 measurements at Arlington",xaxt="n", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(which(NO2_fin[,1]=="1995-02-01"),which(NO2_fin[,1]=="2000-01-01")
             ,which(NO2_fin[,1]=="2005-01-01"),which(NO2_fin[,1]=="2010-01-05")),labels=c(1995,2000,2005,2010), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
Arl_MA=Mov_Avg(Std_Arl, 901, TRUE)
lines(Arl_MA$MA_mean, col="blue", lwd=3,lty=2)
lines(Arl_MA$MA_sd, col="red", lwd=2)

####  Check ACF and PACF
dev.new()
par(mfrow=c(1,2))
acf(Std_Arl)
pacf(Std_Arl)

####  Save data with standardized data added
NO2_fin$Std_Alx=Std_Alx
NO2_fin$Std_Mc=Std_Mc
NO2_fin$Std_Rt=Std_Rt
NO2_fin$Std_Ts=Std_Ts
NO2_fin$Std_Arl=Std_Arl

NO2_fin$Alx_mean=Alx_MA$MA_mean
NO2_fin$Alx_sd=Alx_MA$MA_sd
NO2_fin$Mc_mean=Mc_MA$MA_mean
NO2_fin$Mc_sd=Mc_MA$MA_sd
NO2_fin$Rt_mean=Rt_MA$MA_mean
NO2_fin$Rt_sd=Rt_MA$MA_sd
NO2_fin$Ts_mean=Ts_MA$MA_mean
NO2_fin$Ts_sd=Ts_MA$MA_sd
NO2_fin$Arl_mean=Arl_MA$MA_mean
NO2_fin$Arl_sd=Arl_MA$MA_sd

NewData=NO2_fin
save(NewData,file="NewData.Rdata")
load(file = "NewData.Rdata")