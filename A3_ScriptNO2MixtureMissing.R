####  Load data
library(MASS)
library(plotrix)
source("C:/Linear_prediction/2_TransformedOperations.R")
load(file="C:/Linear_prediction/NewData.Rdata")     #Old data
load(file="C:/Linear_prediction/recentData.Rdata")  #New data

####  Reorder the data
dates <- NewData[,1]
NewData <- NewData[,-1] #mc,rt,ts,arl,alx #train: 1:3442 test: 3443:5163
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
source("C:/Linear_prediction/16_Mix_ECDF_GPD.R")
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

save(NewDat,file="C:/Linear_prediction/NewDat.Rdata")
load(file="C:/Linear_prediction/NewDat.Rdata")

####  Estimate the TPDM
source("C:/Linear_prediction/6_TPDM_Est.R")
TPDM_Out=TPDM_Est(Nrow = 5,X_t = NewDat[,21:25])
TPDM_Out$TPDM_hat
TPDM_Out$TPDM_P_hat

####  Vector b minimizing the tail ratio
TPDM_hat=TPDM_Out$TPDM_hat
b=solve(TPDM_hat[1:4,1:4])%*%TPDM_hat[1:4,5]

####  The best transformed linear predictor
Xhat=Amul(t(b),t(NewDat[,21:24]))
Xhat=as.vector(Xhat)
summary(Xhat)

####  Decompose the TPDM_Pred or the TPDM_Pred_hat
source("C:/Linear_prediction/7_Ang_CPD.R")
AngCPD_Out=Ang_CPD(TPDM_P_hat = TPDM_Out$TPDM_P_hat,m = 70)

save(AngCPD_Out,file="C:/Linear_prediction/AngCPD_NO2_mix_comb.Rdata")
load(file="C:/Linear_prediction/AngCPD_NO2_mix_comb.Rdata")

####  Check the coverage rate
source("C:/Linear_prediction/8_JointRegion.R")
Out1=JointRegion(X_f = NewDat[,25],Xhat = Xhat,Ang = AngCPD_Out$All_ang,Ang_mass = AngCPD_Out$All_mass,tol=0.0001,plot = T,ray = FALSE,cone = T,xr1 = 0,xr2 = 40,yr1 = 0,yr2 = 40)
Out1$Coverage
#[1] 0.98

####  Plot of KDE for angular densities
library(VGAM)
library(ks)
source("C:/Linear_prediction/9_Ang_kde.R")
kde_cpd_out=Ang_kde(Ang_T = AngCPD_Out$All_ang, Ang_mass = AngCPD_Out$All_mass,bw = T,h = 0.3 ,plot_theta = T, plot_w = T,y_lim = 2)

####  Plot conditional intervals
source("C:/Linear_prediction/10_CondDens.R")
a=2
CondDens(z = Out1$Mtx_P[a,],h_w2_CPD = kde_cpd_out$kde_trans_w,plot = T)
##  Inputs: cbind(X_predictant, Xhat), true kde, kde via CPD
##  Outputs: X_(p+1), conditional density, CDF / X_(p+1), conditional density via CPD, CDF via CPD

####  Pick up dates when large NO2 levels exist
source("C:/Linear_prediction/11_CondInterval.R")
Large1=unlist(recentDat[recentDat$Date=="02/04/2019",c(19,20,21,22)])
Large2=unlist(recentDat[recentDat$Date=="02/05/2019",c(19,20,21,22)])
Large3=unlist(recentDat[recentDat$Date=="03/13/2019",c(19,20,21,22)])
Large4=unlist(recentDat[recentDat$Date=="12/23/2019",c(19,20,21,22)])
Large5=unlist(recentDat[recentDat$Date=="01/23/2020",c(19,20,21,22)])

Xhat_e=c(Amul(t(b),Large1),Amul(t(b),Large2),Amul(t(b),Large3),Amul(t(b),Large4),Amul(t(b),Large5))

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

####  Plot all conditional intervals for the top 5%
dev.new()
par(mar=c(6,5.3,1.5,1.5))
Date_x=c("02/04/2019","02/05/2019","03/13/2019","12/23/2019","01/23/2020")
Date_x=as.Date(Date_x,format = "%m/%d/%Y")
plot(1:5,Xhat_e,ylim=c(0,50),col="blue",xlab="",ylab=expression(NO[2]), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, pch=20, cex=1.5, xaxt="n")
axis(1, 1:5, format(Date_x,"%m/%d/%Y"),las=2,cex.axis=1.3)
points(1:5+0.1,Xhat_p,pch=4,cex=1.5,col="red")

for(a in 1:5){
  ##  Cond dens
  CondDen_Out=CondDens(z = Xhat_e[a],h_w2_CPD = kde_cpd_out$kde_trans_w, plot = FALSE)
  ##  Cond den via CPD
  Pred_inter_CPD=CondInterval(z2 = CondDen_Out$z2_CPD,X_f_single = Xhat_e[a] ,cumTraps = CondDen_Out$cumTraps, tol = 0.0001)
  
  points(x = a,y=Pred_inter_CPD$Con_L,pch="-", col="blue",cex=2,lwd=2)
  points(x = a,y=Pred_inter_CPD$Con_U,pch="-", col="blue",cex=2,lwd=2)
  lines(x = c(a,a),c(Pred_inter_CPD$Con_L,Pred_inter_CPD$Con_U),col="blue",cex=2,lwd=2)
  
  points(x = a+0.1,y=Xhat_p[a]-1.96*sqrt(MSPE),pch="-", col="red",cex=2,lwd=2)
  points(x = a+0.1,y=Xhat_p[a]+1.96*sqrt(MSPE),pch="-", col="red",cex=2,lwd=2)
  lines(x = c(a,a)+0.1,c(Xhat_p[a]-1.96*sqrt(MSPE),Xhat_p[a]+1.96*sqrt(MSPE)),col="red",cex=2,lwd=2)
  
  print(c(Pred_inter_CPD$Con_L,Pred_inter_CPD$Con_U))
}

####  Find prediction intervals for the largest 5%
source("C:/Linear_prediction/14_CondIntervalPlot.R")
Cond_Out=CondIntervalPlot(Mtx_P = Out1$Mtx_P,Thres = 0.95,h_w2_CPD = kde_cpd_out$kde_trans_w)
Cond_mtx=Cond_Out$CondIntervalSave

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
z=Cond_mtx[,1]
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
  ##  Cond dens
  CondDen_Out=CondDens(z = Xhat_e[a],h_w2_CPD = kde_cpd_out$kde_trans_w, plot = FALSE)
  ##  Cond den via CPD
  Pred_inter_CPD=CondInterval(z2 = CondDen_Out$z2_CPD,X_f_single = Xhat_e[a], cumTraps = CondDen_Out$cumTraps, tol = 0.0001)
  
  L_orig=InvFt_mix(Pred_inter_CPD$Con_L)
  U_orig=InvFt_mix(Pred_inter_CPD$Con_U)
  
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

