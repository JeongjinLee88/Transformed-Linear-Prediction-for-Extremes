####  Load data
library(MASS)
library(plotrix)
setwd("/home/leej40/Documents/extlinear/output")
load(file="NewData.Rdata")
head(NewData)
dim(NewData)[1]*(2/3) #3442
setwd("/home/leej40/Documents/extlinear/code")
source("TransformedOperations.R")

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
source("Mix_ECDF_GPD.R")
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

save(Train_NO2,file="Train_NO2_Alx.Rdata")
save(Test_NO2,file="Test_NO2_Alx.Rdata")

load(file="Train_NO2_Alx.Rdata")
load(file="Test_NO2_Alx.Rdata")

####  Estimate the TPDM
source("estimateParams.R")
Thres_u=0.95
Est=estimateParams(X = Train_NO2[,21:25], Thres = Thres_u) # Pareto variables
Est$TPDM_hat
Est$TPDM_Phat

####  Vector b minimizing the tail ratio
TPDM_hat=Est$TPDM_hat
bhat=Est$bhat

####  The best transformed linear predictor
Xhat=Amul(t(bhat),t(Test_NO2[,21:24]))
Xhat=as.vector(Xhat)
summary(Xhat)

####  Decompose the TPDM_Pred or the TPDM_Pred_hat
source("CPfactor.R")
CP_no2=CPfactor(Mtx = Est$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 50)
save(CP_no2,file="CP_no2.RData")
load(file="CP_no2.RData")

####  Check the coverage rate
####  Create the 95% joint polar region from a CP-factorization
library(plotrix)  # 'draw.circle'
source("jointRegion.R")
##  Find the 95% joint polar region (Figure 3 (left))
jointOut=jointRegion(Xhat = Xhat, Xf = Test_NO2[,25],
                     Angular = CP_no2$angular, Pmass = CP_no2$pmass, Quan = 0.95,
                     Plot = T, axisLimit = 40, dataPoint = 1)
jointOut$coverage

##  Find the bandwidth via cross-validation
seqbw=seq(0.1,0.7,by=0.05) # a seq of bandwidth
cv_bw=rep(NA,length(seqQuan))
## coverage rates for each target quantile
cv_result=sapply(seqbw,function(bw) crossValidate(Dat = rbind(Train_NO2[,21:25],Test_NO2[,21:25]),Ang =CP_no2$angular,
                                                  pMass = CP_no2$pmass,Thres = Thres_u,
                                                  kfold = 3,
                                                  bandW = bw,Quan = 0.95))
##  Find cv-bandwidths
cv_bw=seqbw[which.min(abs(cv_result-0.95))]


####  Plot of KDE for angular densities
library(VGAM)
library(ks)
source("KDE_w.R")
kde_out=KDE_w(Ang = CP_no2$angular,Pmass = CP_no2$pmass,bw = T,h = 0.5,Plot=T)
kde_out

####  Plot an approximate conditional density with the 95% conditional interval
source("condDensity.R")
conden_out=condDensity(xhatPoint=Xhat_test[471],xp1Point=Test[471,7],
                       kde_h = kde_out, Quan = 0.95,
                       xlim_upper = Xhat_test[471]+0.05, ylim_upper = 100, Plot_h = FALSE)

####  Assess the coverage rate
####  Plot conditional intervals with lines to reproduce Figure 4 (right)
source("coverageRate.R")
target_rate=0.95
XhatXp1 <- cbind(Xhat,Test_NO2[,25])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = target_rate,Plot = T)
coverOut$CoverageRate # 0.965
dev.off()


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

z=coverOut$PI[,1] # lower bound

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
Rad <- jointOut$XY[,1] # Xhat in a test set
W <- jointOut$XY[,1]/Rad
Keep <- Rad > quantile(Rad, .95)  
Test_top5=Test_NO2[Keep,]
Mtx_P_top=XhatXp1_top # two columns (Xhat, Xp+1) given large Xhat
Xhat_o=InvFt_mix(z = Mtx_P_top[,1])*Test_top5[,"Alx_sd"]+Test_top5[,"Alx_mean"]
X_f_o=InvFt_mix(z = Mtx_P_top[,2])*Test_top5[,"Alx_sd"]+Test_top5[,"Alx_mean"]

#Xhat_o=InvFt_mix(z = Mtx_P_top[,1])*Train_top5[,"Alx_sd"]+Train_top5[,"Alx_mean"]
#X_f_o=InvFt_mix(z = Mtx_P_top[,2])*Train_top5[,"Alx_sd"]+Train_top5[,"Alx_mean"]

Mtx_orig=cbind(Xhat_o,X_f_o)
summary(Mtx_orig)

##  Transform the scale of the conditional intervals back to the original scale
L_orig=InvFt_mix(coverOut$PI[,1])
U_orig=InvFt_mix(coverOut$PI[,2])

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
pdf("/home/leej40/Documents/extlinear/no2OriginalScale.pdf",6,6)
par(mar=c(5.1,5.1,2,2))
plot(Mtx_orig[,1],Mtx_orig[,2],xlim=c(0,max(Mtx_orig[,2])*1.3),ylim=c(0,max(Mtx_orig[,2])*1.3),
     main="",xlab=expression(paste("Predicted"," ",NO[2]," ","(ppb)")),ylab=expression(paste("Observed"," ",NO[2]," ","(ppb)")), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.5)

for(a in 1:length(Mtx_orig[,1])){
  points(x = Mtx_orig[a,1],y=L_o[a],pch="-", col="blue",cex=2,lwd=2)
  points(x = Mtx_orig[a,1],y=U_o[a],pch="-", col="blue",cex=2,lwd=2)
}
dev.off()

####  Check the coverage rate with the original scale
sum(L_o <= Mtx_orig[,2] & Mtx_orig[,2] <= U_o)/length(Mtx_orig[,2])
#[1] 0.965    ## Same