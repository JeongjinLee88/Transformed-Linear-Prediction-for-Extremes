####  Load data
library(MASS)
library(plotrix)
source("C:/Linear_prediction/2_TransformedOperations.R")
load(file="C:/Linear_prediction/NewData.Rdata")
head(NewData)
dim(NewData)[1]*(2/3) #3442

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

####  Divide the data set into two sets randomly
set.seed(12345)
R_index=sample(x = seq(1:5163),size = 5163,replace = F)
NewData=NewData[R_index,]

Train_NO2=NewData[1:3442,]
Test_NO2=NewData[3443:5163,]

####  Find a mixture dist of an empirical CDF and a GPD
library(evd)
library(ismev)
library(Hmisc)
source("C:/Linear_prediction/16_Mix_ECDF_GPD.R")
shift=0.9352074
#mrl.plot(Train_NO2[,10])  #Standardized data at Arl
#quantile(Train_NO2[,10],probs = 0.95)
# 1.758

Pa_Train=matrix(0,dim(Train_NO2)[1],5)
Pa_Test=matrix(0,dim(Test_NO2)[1],5)
for(a in 1:5){
  ##  Fit a training set to a GPD
  Gpd_Out=gpd.fit(xdat = Train_NO2[,a+5],threshold = quantile(Train_NO2[,a+5],probs = 0.95))
  scale = Gpd_Out$mle[1]
  shape = Gpd_Out$mle[2]
  rate =  as.numeric(Gpd_Out$rate)
  thold = as.numeric(Gpd_Out$threshold)
  
  ##  Find empirical CDFs under the threshold
  Fn=ECDF(Train_NO2[,a+5])
  CDF=Fn(Train_NO2[,a+5])
  U <- CDF
  sl <- sort.list(Train_NO2[,a+5])
  Train_NO2_sort <- Train_NO2[,a+5][sl]
  U_sort <- U[sl]
  below <- Train_NO2_sort < thold
  CDF_c=cbind(Train_NO2_sort[below],U_sort[below])
  
  ##  Save the results
  Pa_Train[,a]=Mix_ECDF_GPD(Train_NO2[,a+5])
  Pa_Test[,a]=Mix_ECDF_GPD(Test_NO2[,a+5])
}

summary(Pa_Train)
summary(Pa_Test)

Train_NO2=cbind(Train_NO2,Pa_Train)
#Train_NO2=na.omit(Train_NO2)
Test_NO2=cbind(Test_NO2,Pa_Test)
#Test_NO2=na.omit(Test_NO2)

####  Check the summary stat
summary(Train_NO2[,21]) #2/(1+shift)
summary(InvT(x = Train_NO2[,21]))  # zero

save(Train_NO2,file="C:/Linear_prediction/Train_NO2_Alx.Rdata")
save(Test_NO2,file="C:/Linear_prediction/Test_NO2_Alx.Rdata")

load(file="C:/Linear_prediction/Train_NO2_Alx.Rdata")
load(file="C:/Linear_prediction/Test_NO2_Alx.Rdata")

####  Estimate the TPDM
source("C:/Linear_prediction/6_TPDM_Est.R")
TPDM_Out=TPDM_Est(Nrow = 5,X_t = Train_NO2[,21:25]) #Take variables on Pareto scales
TPDM_Out$TPDM_hat
TPDM_Out$TPDM_P_hat

####  Vector b minimizing the tail ratio
TPDM_hat=TPDM_Out$TPDM_hat
b=solve(TPDM_hat[1:4,1:4])%*%TPDM_hat[1:4,5]

####  The best transformed linear predictor
Xhat=Amul(t(b),t(Test_NO2[,21:24]))
Xhat=as.vector(Xhat)
summary(Xhat)

####  Decompose the TPDM_Pred or the TPDM_Pred_hat
source("C:/Linear_prediction/7_Ang_CPD.R")
AngCPD_Out=Ang_CPD(TPDM_P_hat = TPDM_Out$TPDM_P_hat,m = 70)

save(AngCPD_Out,file="C:/Linear_prediction/AngCPD_NO2_mix_Alx.Rdata")
load(file="C:/Linear_prediction/AngCPD_NO2_mix_Alx.Rdata")

####  Check the coverage rate
source("C:/Linear_prediction/8_JointRegion.R")
Out1=JointRegion(X_f = Test_NO2[,25],Xhat = Xhat,Ang = AngCPD_Out$All_ang,Ang_mass = AngCPD_Out$All_mass,tol=0.0001,plot = T,ray = FALSE,cone = T,xr1 = 0,xr2 = 40,yr1 = 0,yr2 = 40)
Out1$Coverage
#[1] 0.988

####  Plot of KDE for angular densities
library(VGAM)
library(ks)
source("C:/Linear_prediction/9_Ang_kde.R")
kde_cpd_out=Ang_kde(Ang_T = AngCPD_Out$All_ang, Ang_mass = AngCPD_Out$All_mass,bw = T,h = 0.3 ,plot_theta = T, plot_w = T,y_lim = 2)

####  Conditional densities
source("C:/Linear_prediction/10_CondDens.R")
a=2
CondDen_Out=CondDens(z = Out1$Mtx_P[a,],h_w2_CPD = kde_cpd_out$kde_trans_w, plot = T)
##  Inputs: cbind(X_predictant, Xhat), true kde, kde via CPD
##  Outputs: X_(p+1), conditional density, CDF / X_(p+1), conditional density via CPD, CDF via CPD

####  Coverage rate for the conditional density
source("C:/Linear_prediction/11_CondInterval.R")
Pred_inter=CondInterval(z2 = CondDen_Out$z2_CPD, X_f_single = Out1$Mtx_P[a,2], cumTraps = CondDen_Out$cumTraps, tol = 0.0001)

####  Check the coverage rate in a 'test' set
####  Conditional densities for the largest 5% determined by Xhat
####  Reproduce Figure 5 in air pollution application
source("C:/Linear_prediction/13_AssessCoverage.R")
Assess_out=AssessCoverage(Mtx_P = Out1$Mtx_P,h_w2_CPD = kde_cpd_out$kde_trans_w, Thres = 0.95,Thres_Xhat = T)
#[1] 0.965

####  Check the coverage rate for the largest 5% for which ||X1,X2,X3,X4,Xhat|| is large
Mtx_obs=cbind(Test_NO2[,21:24],Xhat) #X1,X2,X3,X4,Xhat in a test set
Rad1 <- sqrt(apply(Mtx_obs^2,1,sum))
Keep1 <- Rad1 > quantile(Rad1, .95)  
Mtx_obs[Keep1,]

Rad <- Mtx_obs[,5] #Xhat in a test set
Keep <- Rad > quantile(Rad, .95)  
Mtx_obs[Keep,]
dim(Mtx_obs[Keep1&Keep,]) #66 out of 86

Rad <- Mtx_obs[Keep1,5] #Xhat for which ||X1,X2,X3,X4,Xhat|| is large in a test set
W <- Test_NO2[Keep1,25]/Rad
Mtx_P_top=cbind(Mtx_obs[Keep1,5],Test_NO2[Keep1,25])
AssessCoverage(Mtx_P = Mtx_P_top,h_w2_CPD = kde_cpd_out$kde_trans_w, Thres = 0.95,Thres_Xhat = T)
#[1] 0.9534

####  Plot all conditional intervals for the top 5% data
source("C:/Linear_prediction/14_CondIntervalPlot.R")
Cond_out=CondIntervalPlot(Mtx_P = Out1$Mtx_P,Thres = 0.95,h_w2_CPD = kde_cpd_out$kde_trans_w)
Cond_mtx=Cond_out$CondIntervalSave # (LB,UB,Xhat,X_f)

####  Plot conditional intervals with lines
dev.new()
par(mar=c(5.1,5.1,2,2))
plot(Cond_mtx[,3],Cond_mtx[,4],xlim=c(0,55),ylim=c(0,55),
     main="",xlab=expression(paste("Predicted"," ",NO[2]," ","(Pareto Scale)")),ylab=expression(paste("Observed"," ",NO[2]," ","(Pareto Scale)")), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.5)
sl=sort.list(Cond_mtx[,3])
lines(Cond_mtx[sl,3],Cond_mtx[sl,1], lty = 2, cex=2, lwd=2, col="blue")
lines(Cond_mtx[sl,3],Cond_mtx[sl,2], lty = 2, cex=2, lwd=2, col="blue")

################################
####  Transform back to original scale
################################

####  Back to the original scale
Gpd_Out=gpd.fit(xdat = Train_NO2[,"Std_Alx"],threshold = quantile(Train_NO2[,"Std_Alx"],probs = 0.95))
scale = Gpd_Out$mle[1]
shape = Gpd_Out$mle[2]
rate =  as.numeric(Gpd_Out$rate)
thold = as.numeric(Gpd_Out$threshold)
tholdInv <- qgpd(1 - rate, loc = 1, scale = 0.5, shape = 0.5)

#overlap <- 1
origData <- Train_NO2[,"Std_Alx"]
Fn=ECDF(Train_NO2[,"Std_Alx"])
CDF=Fn(Train_NO2[,"Std_Alx"])
uData <- CDF
sl <- sort.list(origData)
sortOrigData <- origData[sl]
sortUData <- uData[sl]
below <- sortOrigData < thold
#below <- sortOrigData < thold + overlap
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
    
    #u <- 1-(z[below]+shift)^(-2)
    #s <- seq(1, length(z[below]))
    #positive <- s[u > 0 & u < 1]
    #negative <- s[u <= 0 & !is.na(u)]
    #x[positive] <- approxExtrap(x = CDF_c[,2],y = CDF_c[,1],xout = u[positive])$y
    #u[u<0]=0.1
    #x[below] <- approxExtrap(x = CDF_c[,2],y = CDF_c[,1],xout = u)$y
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

####  Choose the top 5%
Rad <- Out1$Mtx_P[,1] # Xhat in a test set
W <- Out1$Mtx_P/Rad
Keep <- Rad > quantile(Rad, .95)  
Test_top5=Test_NO2[Keep,]
#Train_top5=Train_NO2[Keep,]
Mtx_P_top=Assess_out$Mtx_P_top
Xhat_o=InvFt_mix(z = Mtx_P_top[,1])*Test_top5[,"Alx_sd"]+Test_top5[,"Alx_mean"]
X_f_o=InvFt_mix(z = Mtx_P_top[,2])*Test_top5[,"Alx_sd"]+Test_top5[,"Alx_mean"]

#Xhat_o=InvFt_mix(z = Mtx_P_top[,1])*Train_top5[,"Alx_sd"]+Train_top5[,"Alx_mean"]
#X_f_o=InvFt_mix(z = Mtx_P_top[,2])*Train_top5[,"Alx_sd"]+Train_top5[,"Alx_mean"]

Mtx_orig=cbind(Xhat_o,X_f_o)
summary(Mtx_orig)

##  Transform the scale of the conditional intervals back to the original scale
L_orig=InvFt_mix(Cond_mtx[,1])
U_orig=InvFt_mix(Cond_mtx[,2])

L_o=L_orig*Test_top5[,"Alx_sd"]+Test_top5[,"Alx_mean"]
U_o=U_orig*Test_top5[,"Alx_sd"]+Test_top5[,"Alx_mean"]

Int_orig_linear=cbind(L_o,U_o)
#L_o=L_orig*Train_top5[,"Alx_sd"]+Train_top5[,"Alx_mean"]
#U_o=U_orig*Train_top5[,"Alx_sd"]+Train_top5[,"Alx_mean"]

summary(L_o)
summary(U_o)

####  All conditional intervals for the top 5% on the original scale
####  Reproduce Figure 5
dev.new()
par(mar=c(5.1,5.1,2,2))
#plot(Mtx_orig[,1],Mtx_orig[,2],xlim=c(0,max(Mtx_orig[,2])*1.5),ylim=c(0,max(Mtx_orig[,2])*1.5),
#     main=expression(paste("The scatter plot of"," ",X[p+1]," ","and"," ",hat(X)[p+1]," ","with conditional intervals")),xlab=expression(hat(X)[p+1]),ylab=expression(X[p+1]), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,pch=20,cex=1)
plot(Mtx_orig[,1],Mtx_orig[,2],xlim=c(0,max(Mtx_orig[,2])*1.3),ylim=c(0,max(Mtx_orig[,2])*1.3),
     main="",xlab=expression(paste("Predicted"," ",NO[2]," ","(ppb)")),ylab=expression(paste("Observed"," ",NO[2]," ","(ppb)")), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.5)

for(a in 1:length(Mtx_orig[,1])){
  points(x = Mtx_orig[a,1],y=L_o[a],pch="-", col="blue",cex=2,lwd=2)
  points(x = Mtx_orig[a,1],y=U_o[a],pch="-", col="blue",cex=2,lwd=2)
}

####  Check the coverage rate with the original scale
sum(L_o <= Mtx_orig[,2] & Mtx_orig[,2] <= U_o)/length(Mtx_orig[,2])
#[1] 0.965    ## Same