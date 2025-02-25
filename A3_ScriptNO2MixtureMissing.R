setwd("/home/leej40/Documents/extlinear/output")
####  Load data
library(MASS)
library(plotrix)
load(file="NewData.Rdata")     
# Alx, Mc, Rt, Ts, Arl
# 02/01/1995 - 03/20/2014
load(file="recentData.RData")  
# Mc, Rt, Ts, Arl / no Alx
# 06/07/2016 - 12/31/2020
setwd("/home/leej40/Documents/extlinear/code")
source("TransformedOperations.R")

####  Reorder the data
#mc,rt,ts,arl,alx #train: 1:3442 test: 3443:5163
dates <- NewData[,1]
NewData <- NewData[,-1] 
names(NewData[,c(2,3,4,5,1,7,8,9,10,6,11:20)])
NewData=NewData[,c(2,3,4,5,1,7,8,9,10,6,11:20)]
dim(NewData)
##  Columns: 1-5 original data, 6-10 standardized data, 11-20 MA

####  Define an empirical CDF
ECDF<-function (x){
  x <- sort(x)
  n <- length(x)
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/(n+1), 
                    method = "constant", yleft = 0, yright = 1, f = 0, 
                    ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}

####  The old data are considered as a training data
####  and the new data as a test set.
set.seed(12345)
R_index=sample(x = seq(1:5163),size = 5163,replace = F)
NewData=NewData[R_index,]

shift=0.9352074
Pa_train_nonpa=matrix(0,nrow = dim(NewData)[1],ncol = 5)
for(i in 1:5){
  Fn=ECDF(NewData[1:dim(NewData)[1],i+5])
  Pa_train_nonpa[,i]=(1-Fn(NewData[1:dim(NewData)[1],i+5]))^(-1/2)-shift
}

N_dist=matrix(0,nrow=dim(NewData)[1],ncol=5)
for(i in 1:5){
  Fn=ECDF(NewData[,i+5])
  N_dist[,i]=qnorm(p = Fn(NewData[,i+5]),mean = 0,sd = 1)
}

N_recent=matrix(0,nrow=dim(recentDat)[1],ncol=4)
for(i in 1:4){
  Fn=ECDF(recentDat[,i+5])
  N_recent[,i]=qnorm(p = Fn(recentDat[,i+5]),mean = 0,sd = 1)
}

####  A mixture dist
library(evd)
library(ismev)
library(Hmisc)
source("Mix_ECDF_GPD.R")
shift=0.9352074

Pa_old=matrix(0,dim(NewData)[1],5)
Pa_new=matrix(0,dim(recentDat)[1],5)
for(a in 1:5){
  Gpd_Out=gpd.fit(xdat = NewData[,a+5],threshold = quantile(NewData[,a+5],probs = 0.95))
  scale = Gpd_Out$mle[1]
  shape = Gpd_Out$mle[2]
  rate =  as.numeric(Gpd_Out$rate)
  thold = as.numeric(Gpd_Out$threshold)
  
  origData <- NewData[,a+5]
  Fn=ECDF(origData)
  CDF=Fn(origData)
  uData <- CDF
  sl <- sort.list(origData)
  sortOrigData <- origData[sl]
  sortUData <- uData[sl]
  below <- sortOrigData < thold
  CDF_c=cbind(sortOrigData[below],sortUData[below])
  
  ##  Save the results
  Pa_old[,a]=Mix_ECDF_GPD(NewData[,a+5])
  Pa_new[,a]=Mix_ECDF_GPD(recentDat[,a+5])
}

##  There are no obs at Alx in "Pa_new"
Pa_new=Pa_new[,1:4]
NewDat=cbind(NewData,Pa_old)
recentDat=cbind(recentDat,Pa_new)

save(NewDat,file="NewDat.Rdata")
load(file="/home/leej40/Documents/extlinear/output/NewDat.Rdata")

####  Estimate the TPDM
source("estimateParams.R")
Thres_u=0.95
Est=estimateParams(X = NewDat[,21:25], Thres = Thres_u) # Pareto variables
Est$TPDM_hat
Est$TPDM_Phat

bhat=Est$bhat
####  The best transformed linear predictor
Xhat=Amul(t(bhat),t(NewDat[,21:24]))
Xhat=as.vector(Xhat)
summary(Xhat)

####  Decompose the TPDM_Pred or the TPDM_Pred_hat
source("CPfactor.R")
CP_no2_comb=CPfactor(Mtx = Est$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 50)
save(CP_no2_comb,file="CP_no2_comb.RData")
load(file="CP_no2_comb.RData")

####  Check the coverage rate
####  Create the 95% joint polar region from a CP-factorization
library(plotrix)  # 'draw.circle'
source("jointRegion.R")
##  Find the 95% joint polar region (Figure 3 (left))
jointOut=jointRegion(Xhat = Xhat, Xf = NewDat[,25],
                     Angular = CP_no2_comb$angular, Pmass = CP_no2_comb$pmass, Quan = 0.95,
                     Plot = T, axisLimit = 40, dataPoint = 1)
jointOut$coverage

##  Find the bandwidth via cross-validation
seqbw=seq(0.1,0.7,by=0.05) # a seq of bandwidth
cv_bw=rep(NA,length(seqQuan))
## coverage rates for the target quantile
cv_result=sapply(seqbw,function(bw) crossValidate(Dat = NewDat[,21:25],Ang =CP_no2_comb$angular,
                                                  pMass = CP_no2_comb$pmass,Thres = Thres_u,
                                                  kfold = 3,
                                                  bandW = bw,Quan = 0.95))
##  Find cv-bandwidths
cv_bw=seqbw[which.min(abs(cv_result-0.95))]

####  Plot of KDE for angular densities
library(VGAM)
library(ks)
source("KDE_w.R")
kde_out=KDE_w(Ang = CP_no2_comb$angular,Pmass = CP_no2_comb$pmass,bw = T,h = cv_bw, Plot=T)
kde_out

####  Plot an approximate conditional density with the 95% conditional interval
source("condDensity.R")
largeObs=cbind(Xhat,NewDat[,25])
conden_out=condDensity(xhatPoint=largeObs[1],xp1Point=largeObs[2],
                       kde_h = kde_out, Quan = 0.95,
                       xlim_upper = largeObs[1]+0.07, ylim_upper = 70, Plot_h = FALSE)


####  Pick up dates when large NO2 levels exist
Large1=unlist(recentDat[recentDat$Date=="02/04/2019",c(19,20,21,22)])
Large2=unlist(recentDat[recentDat$Date=="02/05/2019",c(19,20,21,22)])
Large3=unlist(recentDat[recentDat$Date=="03/13/2019",c(19,20,21,22)])
Large4=unlist(recentDat[recentDat$Date=="12/23/2019",c(19,20,21,22)])
Large5=unlist(recentDat[recentDat$Date=="01/23/2020",c(19,20,21,22)])

Xhat_e=c(Amul(t(bhat),Large1),Amul(t(bhat),Large2),Amul(t(bhat),Large3),Amul(t(bhat),Large4),Amul(t(bhat),Large5))

####  Gaussian cases
Cov=cov(N_dist)
b_n=solve(Cov[1:4,1:4])%*%Cov[1:4,dim(N_dist)[2]]
Xhat_n=t(t(b_n)%*%t(N_recent))
MSPE=Cov[5,5]-t(Cov[1,2:5])%*%solve(Cov[1:4,1:4])%*%Cov[1,2:5]
MSPE=as.numeric(MSPE)
N_recent_date=cbind.data.frame(recentDat[,1],N_recent)
Large1=N_recent_date[N_recent_date[,1]=="02/04/2019",2:5]
Large2=N_recent_date[N_recent_date[,1]=="02/05/2019",2:5]
Large3=N_recent_date[N_recent_date[,1]=="03/13/2019",2:5]
Large4=N_recent_date[N_recent_date[,1]=="12/23/2019",2:5]
Large5=N_recent_date[N_recent_date[,1]=="01/23/2020",2:5]

Xhat_p=c(t(b_n)%*%t(Large1),t(b_n)%*%t(Large2),t(b_n)%*%t(Large3),t(b_n)%*%t(Large4),t(b_n)%*%t(Large5))
#S2_pred1=MSPE*(1+t(Large1)%*%solve(t(N_dist[,1:4])%*%N_dist[,1:4])%*%Large1)

##  coverage rates for all large Xhat based on Gaussian dist
##  1. Standardized data
Xhat_n=t(t(b_n)%*%t(NewData[,6:9]))
Comb_n=cbind(Xhat_n,NewData[,10])
##  2. Square-root
Xhat_n=t(t(b_n)%*%t(sqrt(NewData[,1:4])))
Comb_n=cbind(Xhat_n,sqrt(NewData[,5]))
##  3. Original scale
Xhat_n=t(t(b_n)%*%t(NewData[,1:4]))
Comb_n=cbind(Xhat_n,NewData[,5])
Q_n=quantile(x = Xhat_n,probs = 0.95)
Top_n=Comb_n[Comb_n[,1]>Q_n,]

UB=Top_n[,1]+1.96*sqrt(MSPE)
LB=Top_n[,1]-1.96*sqrt(MSPE)

Count_n=Top_n[,2]>LB & Top_n[,2]<UB
sum(Count_n)/length(Count_n)

####  Plot conditional intervals for five dates on Pareto scale
dev.new()
par(mar=c(6,5.3,1.5,1.5))
Date_x=c("02/04/2019","02/05/2019","03/13/2019","12/23/2019","01/23/2020")
Date_x=as.Date(Date_x,format = "%m/%d/%Y")
plot(1:5,Xhat_e,ylim=c(0,50),col="blue",xlab="",ylab=expression(NO[2]), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, pch=20, cex=1.5, xaxt="n")
axis(1, 1:5, format(Date_x,"%m/%d/%Y"),las=2,cex.axis=1.3)
points(1:5+0.1,Xhat_p,pch=4,cex=1.5,col="red")

for(a in 1:5){
  CondDen_Out=condDensity(xhatPoint=Xhat_e[a],xp1Point = NA,
                kde_h = kde_out, Quan = 0.95,
                xlim_upper = largeObs[1]+0.07, ylim_upper = 70,Plot_v = F, Plot_h = FALSE)
  points(x = a,y=CondDen_Out$condQuan[1],pch="-", col="blue",cex=2,lwd=2)
  points(x = a,y=CondDen_Out$condQuan[2],pch="-", col="blue",cex=2,lwd=2)
  lines(x = c(a,a),CondDen_Out$condQuan,col="blue",cex=2,lwd=2)
  
  points(x = a+0.1,y=Xhat_p[a]-1.96*sqrt(MSPE),pch="-", col="red",cex=2,lwd=2)
  points(x = a+0.1,y=Xhat_p[a]+1.96*sqrt(MSPE),pch="-", col="red",cex=2,lwd=2)
  lines(x = c(a,a)+0.1,c(Xhat_p[a]-1.96*sqrt(MSPE),Xhat_p[a]+1.96*sqrt(MSPE)),col="red",cex=2,lwd=2)
}

####  Find prediction intervals for the largest 5%
source("coverageRate.R")
target_rate=0.95
XhatXp1 <- cbind(Xhat,NewDat[,25])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = target_rate,Plot = T)
coverOut$CoverageRate

####  All conditional intervals on the original scale
Gpd_Out=gpd.fit(xdat = NewDat[,"Std_Alx"],threshold = quantile(NewDat[,"Std_Alx"],probs = 0.95))
scale = Gpd_Out$mle[1]
shape = Gpd_Out$mle[2]
rate =  as.numeric(Gpd_Out$rate)
thold = as.numeric(Gpd_Out$threshold)
tholdInv <- qgpd(1 - rate, loc = 1, scale = 0.5, shape = 0.5)

origData <- NewDat[,"Std_Alx"]
Fn=ECDF(NewDat[,"Std_Alx"])
CDF=Fn(NewDat[,"Std_Alx"])
uData <- CDF
sl <- sort.list(origData)
sortOrigData <- origData[sl]
sortUData <- uData[sl]
below <- sortOrigData < thold
CDF_c=cbind(sortOrigData[below],sortUData[below])
#Fn=ECDF(sortOrigData[below])
#CDF=Fn(sortOrigData[below])
#CDF_c=cbind(sortOrigData[below],CDF)
#CDF_c=unique(CDF_c)
#xHold=seq(min(sortOrigData[below]), max(sortOrigData[below]), 0.01)
#uHold=approx(x = CDF_c[,1],y = CDF_c[,2],xout = xHold)$y
#uHold=approxExtrap(x = CDF_c[,1],y = CDF_c[,2],xout = xHold)$y
z=coverOut$PI[,1]
InvFt_mix <- function(z)
{
  x <- numeric(length(z))
  below <- z <= tholdInv
  
  ##  below thold
  if(length(z[below] > 0))
  {
    xBelow <- x[below]
    u <- 1-(z[below]+shift)^(-2)
    IR <- u < 1 & u > 0 & !is.na(u)
    xBelow[IR] <- approxExtrap(x = CDF_c[,2],y = CDF_c[,1],xout = u[IR])$y
    xBelow[!IR] <- NA
    x[below]<-xBelow
    
    ##u <- 1-(z[below]+shift)^(-2)
    #s <- seq(1, length(z[below]))
    #positive <- s[u > 0 & u < 1]
    #negative <- s[u <= 0 & !is.na(u)]
    #x[positive] <- approxExtrap(x = CDF_c[,2],y = CDF_c[,1],xout = u[positive])$y
    ##u[u<0]=0.1
    ##x[below] <- approxExtrap(x = CDF_c[,2],y = CDF_c[,1],xout = u)$y
    #x[below] <- approxExtrap(uHold, xHold, u)$y
  }
  
  ##  above thold, use gpd inv
  if(length(z[!below] > 0))
  {
    intoFun <- 1 - 1/rate*(1 - pgev(z[!below]+shift,1,0.5,0.5))
    xAbove <- x[!below]
    xAbove[intoFun < 1] <- qgpd(intoFun[intoFun < 1], thold,
                                scale, shape)
    xAbove[intoFun == 1] <- Inf
    x[!below] <- xAbove
  }
  
  return(x)
}

integrand=function(x){(1/(2*sqrt(pi)))*exp(-(1/2)*(x^2))}

##  Convert Gaussian scale to original scale.
InvFt_N=function(z){
  origData <- NewData[,"Std_Alx"]
  Fn=ECDF(NewData[,"Std_Alx"])
  CDF=Fn(NewData[,"Std_Alx"])
  uData <- CDF
  sl <- sort.list(origData)
  sortOrigData <- origData[sl]
  sortUData <- uData[sl]
  CDF_c=cbind(sortOrigData,sortUData)
  CDF_c=unique(CDF_c)
  #u=Fn(z)
  u=integrate(integrand,-Inf,z)
  u=u$value
  IR <- u < 1 & u > 0 & !is.na(u)
  Xhat_N=qnorm(u,0,1)
  #Xhat_N=approx(x = CDF_c[,2],y = CDF_c[,1],xout = u[IR])$y
  Xhat_N[!IR]=Inf
  
  return(Xhat_N)
}


####  Choose large NO2 levels on specific dates
head(recentDat)
Date_p=c("02/04/2019", "02/05/2019", "03/13/2019", "12/23/2019", "01/23/2020")
Xhat_o=rep(NA,5)
for(i in 1:5){
  Xhat_o[i]=InvFt_mix(z = Xhat_e[i])*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  print(InvFt_mix(z = Xhat_e[i]))
}
Xhat_o

Xhat_o_N=rep(NA,5)
for(i in 1:5){
  #Xhat_o_N[i]=InvFt_N(z=Xhat_p[i])*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  Xhat_o_N[i]=Xhat_p[i]*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  #print(InvFt_N(z=Xhat_p[i]))  
}
Xhat_o_N

####  Plot conditional intervals for particular cases
####  Reproduce Figure 5 (right)
dev.new()
pdf("/home/leej40/Documents/extlinear/no2RecentData.pdf",6,6)
par(mar=c(6,5.3,1.5,1.5))
Date_x=c("02/04/2019","02/05/2019","03/13/2019","12/23/2019","01/23/2020")
Date_x=as.Date(Date_x,format = "%m/%d/%Y")
plot(1:5,Xhat_o,ylim=c(0,70),col="blue",xlab="",ylab=expression(paste(NO[2]," ","(ppb)")), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, pch=20, cex=1.5, xaxt="n")
#axis(1, 1:5, format(Date_x,"%m/%d/%Y"),las=2,cex.axis=1.3)
axis(1,labels=FALSE,las=2,cex.axis=1.3)
text(x = 1:length(Date_x)+0.3,
     y = par("usr")[3] - 0.45,
     labels = Date_x,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,adj = 1.1,
     cex = 1.3)
points(1:5+0.1,Xhat_o_N,pch=4,cex=1.5,col="red")

for(a in 1:5){
  CondDen_Out=condDensity(xhatPoint=Xhat_e[a],xp1Point = NA,
                          kde_h = kde_out, Quan = 0.95,
                          xlim_upper = largeObs[1]+0.07, ylim_upper = 70,Plot_v = F, Plot_h = FALSE)
  
  L_orig=InvFt_mix(CondDen_Out$condQuan[1])
  U_orig=InvFt_mix(CondDen_Out$condQuan[2])
  
  L_o=L_orig*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  U_o=U_orig*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  
  points(x = a,y=L_o,pch="-", col="blue",cex=2,lwd=2)
  points(x = a,y=U_o,pch="-", col="blue",cex=2,lwd=2)
  lines(x = c(a,a),c(L_o,U_o),col="blue",cex=2,lwd=2)
  
  #L_orig_N=InvFt_N(z=Xhat_p[a]-1.96*sqrt(MSPE))
  #U_orig_N=InvFt_N(z=Xhat_p[a]+1.96*sqrt(MSPE)) 
  
  L_orig_N=Xhat_p[a]-1.96*sqrt(MSPE)
  U_orig_N=Xhat_p[a]+1.96*sqrt(MSPE)
  
  L_o_N=L_orig_N*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  U_o_N=U_orig_N*recentDat[recentDat$Date==Date_p[i],"Arl_sd"]+recentDat[recentDat$Date==Date_p[i],"Arl_mean"]
  
  points(x = a+0.1,y=L_o_N,pch="-", col="red",cex=2,lwd=2)
  points(x = a+0.1,y=U_o_N,pch="-", col="red",cex=2,lwd=2)
  lines(x = c(a,a)+0.1,c(L_o_N,U_o_N),col="red",cex=2,lwd=2,lty=2)
  
}
dev.off()
