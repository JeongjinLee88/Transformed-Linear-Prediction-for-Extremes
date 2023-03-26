####  Load portfolios data
source("C:/Linear_prediction/2_TransformedOperations.R")
Pf <- read.csv("C:/Portfolios/30_Industry_Portfolios_Daily.csv")
length(Pf[,1]) #24854
Pf[Pf==-99.99]<-NA
Pf=na.omit(Pf)

####  Select data for 1950-2020
Pf <- Pf[Pf[,1] > 19500000,]
dates <- Pf[,1]
Pf <- Pf[, -1]
d <- dim(Pf)[2]
n <- dim(Pf)[1]  #17911

####  Reorder variables
NewPf=cbind(Pf$Food,Pf$Smoke,Pf$Games,Pf$Books,Pf$Hshld,Pf$Clths,Pf$Hlth,Pf$Chems,Pf$Txtls
            ,Pf$Cnstr,Pf$Steel,Pf$FabPr,Pf$ElcEq,Pf$Autos,Pf$Carry,Pf$Mines,Pf$Oil,Pf$Util,Pf$Telcm
            ,Pf$Servs,Pf$BusEq,Pf$Trans,Pf$Whlsl,Pf$Rtail,Pf$Meals,Pf$Fin,Pf$Other,Pf$Coal,Pf$Beer,Pf$Paper)
ColName=c(colnames(Pf[,-c(18,2,24)]),colnames(Pf[,c(18,2,24)]))
colnames(NewPf)=ColName

####  Negate data, set negatives (gains) to zero, and bound away from zero
negPf <- -NewPf
negPf[negPf < 0] <- 0
negPf <- Trans(data.matrix(negPf))
set.seed(1234)
R_index=sample(x = seq(1:dim(negPf)[1]),size = dim(negPf)[1],replace = F)
negPf=negPf[R_index,]
dim(negPf)[1]*(2/3)  #11940

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

##  Find the shift with the combined data
shift=0.9352074
Pa=matrix(0,nrow=dim(negPf)[1],ncol=dim(negPf)[2])
for(i in 1:dim(negPf)[2]){
  Fn=ECDF(negPf[,i])
  Pa[,i]=(1-Fn(negPf[,i]))^(-1/2)-shift
}
summary(Pa)

N_dist=negPf
N_dist=matrix(0,nrow=dim(negPf)[1],ncol=dim(negPf)[2])
for(i in 1:dim(negPf)[2]){
  Fn=ECDF(negPf[,i])
  N_dist[,i]=qnorm(p = Fn(negPf[,i]),mean = 0,sd = 1)
}

Pa_train=matrix(0,nrow = 11940,ncol = dim(negPf)[2])
for(i in 1:dim(negPf)[2]){
  Fn=ECDF(negPf[1:11940,i])
  Pa_train[,i]=(1-Fn(negPf[1:11940,i]))^(-1/2)-shift
}
Pa_test=matrix(0,nrow = 17911-11940,ncol = dim(negPf)[2])
for(i in 1:dim(negPf)[2]){
  Fn=ECDF(negPf[1:11940,i])
  CDF=Fn(negPf[1:11940,i])
  CDF_c=cbind(negPf[1:11940,i],CDF)
  sl=sort.list(negPf[1:11940,i])
  CDF_sort=CDF_c[sl,]
  CDF_sort=unique(CDF_sort)
  Fn1=approx(x = CDF_sort[,1],y = CDF_sort[,2],xout = negPf[11941:17911,i])
  Pa_test[,i]=(1-Fn1$y)^(-1/2)-shift
}
summary(Pa_test)
Pa_test=na.omit(Pa_test)

negPf_train=Pa[1:11940,]
negPf_test=Pa[11941:17911,]

N_train=N_dist[1:11940,]
N_test=N_dist[11941:17911,]

save(negPf_train,file="C:/Linear_prediction/negPf_train.Rdata")
save(negPf_test,file="C:/Linear_prediction/negPf_test.Rdata")
save(N_train,file="C:/Linear_prediction/N_train.Rdata")
save(N_test,file="C:/Linear_prediction/N_test.Rdata")

load(file="C:/Linear_prediction/negPf_train.Rdata")
load(file="C:/Linear_prediction/negPf_test.Rdata")

ColName_coal=ColName[c(1:27,29,30,28)]
ColName_beer=ColName[c(1:28,30,29)]
ColName_paper=ColName
####  Three different data sets
Dat_coal_tr=negPf_train[,c(1:27,29,30,28)]
Dat_beer_tr=negPf_train[,c(1:28,30,29)]
Dat_paper_tr=negPf_train

Dat_coal_te=negPf_test[,c(1:27,29,30,28)]
Dat_beer_te=negPf_test[,c(1:28,30,29)]
Dat_paper_te=negPf_test

N_coal_tr=N_train[,c(1:27,29,30,28)]
N_beer_tr=N_train[,c(1:28,30,29)]
N_paper_tr=N_train

N_coal_te=N_test[,c(1:27,29,30,28)]
N_beer_te=N_test[,c(1:28,30,29)]
N_paper_te=N_test

Cov_coal=cov(N_coal_tr)
Cov_beer=cov(N_beer_tr)
Cov_paper=cov(N_paper_tr)

b_coal=solve(Cov_coal[1:29,1:29])%*%Cov_coal[1:29,dim(N_coal_tr)[2]]
b_beer=solve(Cov_beer[1:29,1:29])%*%Cov_beer[1:29,dim(N_coal_tr)[2]]
b_paper=solve(Cov_paper[1:29,1:29])%*%Cov_paper[1:29,dim(N_coal_tr)[2]]

Xhat_n_coal=t(t(b_coal)%*%t(N_coal_te[,1:29]))
Xhat_n_coal=as.vector(Xhat_n_coal)
Xhat_n_beer=t(t(b_beer)%*%t(N_beer_te[,1:29]))
Xhat_n_beer=as.vector(Xhat_n_beer)
Xhat_n_paper=t(t(b_paper)%*%t(N_paper_te[,1:29]))
Xhat_n_paper=as.vector(Xhat_n_paper)

#MSPE_coal=Cov_coal[30,30]-t(Cov_coal[1,2:29])%*%solve(Cov_coal[1:29,1:29])%*%Cov_coal[1,2:29]
MSPE_coal=Cov_coal[30,30]-Cov_coal[30,1:29]%*%solve(Cov_coal[1:29,1:29])%*%Cov_coal[1:29,30]
MSPE_coal=as.numeric(MSPE_coal)
#MSPE_beer=Cov_beer[30,30]-t(Cov_beer[1,2:29])%*%solve(Cov_beer[1:29,1:29])%*%Cov_beer[1,2:29]
MSPE_beer=Cov_beer[30,30]-Cov_beer[30,1:29]%*%solve(Cov_beer[1:29,1:29])%*%Cov_beer[1:29,30]
MSPE_beer=as.numeric(MSPE_beer)
#MSPE_paper=Cov_paper[30,30]-t(Cov_paper[1,2:29])%*%solve(Cov_paper[1:29,1:29])%*%Cov_paper[1,2:29]
MSPE_paper=Cov_paper[30,30]-Cov_paper[30,1:29]%*%solve(Cov_paper[1:29,1:29])%*%Cov_paper[1:29,30]
MSPE_paper=as.numeric(MSPE_paper)

##  Coal
Comb_n=cbind(Xhat_n_coal,N_coal_te[,30])
Q_n=quantile(x = Xhat_n_coal,probs = 0.95)
Top_n=Comb_n[Comb_n[,1]>Q_n,]

UB=Top_n[,1]+1.96*sqrt(MSPE_coal)
LB=Top_n[,1]-1.96*sqrt(MSPE_coal)

Count_n=Top_n[,2]>LB & Top_n[,2]<UB
sum(Count_n)/length(Count_n)
#0.676
##  Beer
Comb_n=cbind(Xhat_n_beer,N_beer_te[,28])
Q_n=quantile(x = Xhat_n_beer,probs = 0.95)
Top_n=Comb_n[Comb_n[,1]>Q_n,]

UB=Top_n[,1]+1.96*sqrt(MSPE_beer)
LB=Top_n[,1]-1.96*sqrt(MSPE_beer)

Count_n=Top_n[,2]>LB & Top_n[,2]<UB
sum(Count_n)/length(Count_n)
#0.341
##  Paper
Comb_n=cbind(Xhat_n_paper,N_paper_te[,28])
Q_n=quantile(x = Xhat_n_paper,probs = 0.95)
Top_n=Comb_n[Comb_n[,1]>Q_n,]

UB=Top_n[,1]+1.96*sqrt(MSPE_paper)
LB=Top_n[,1]-1.96*sqrt(MSPE_paper)

Count_n=Top_n[,2]>LB & Top_n[,2]<UB
sum(Count_n)/length(Count_n)
#0.261
####  Estimate the TPDM
source("C:/Linear_prediction/6_TPDM_Est.R")
TPDM_Out=TPDM_Est(Nrow = dim(Dat_coal_tr)[2],X_t = Dat_coal_tr)
TPDM_Out=TPDM_Est(Nrow = dim(Dat_coal_tr)[2],X_t = Dat_beer_tr)
TPDM_Out=TPDM_Est(Nrow = dim(Dat_coal_tr)[2],X_t = Dat_paper_tr)
TPDM_Out$TPDM_hat
TPDM_Out$TPDM_P_hat

####  Find the optimized Vector b
TPDM_hat=TPDM_Out$TPDM_hat
b=solve(TPDM_hat[1:dim(Dat_coal_tr)[2]-1,1:dim(Dat_coal_tr)[2]-1])%*%TPDM_hat[1:dim(Dat_coal_tr)[2]-1,dim(Dat_coal_tr)[2]]
which(b==max(b))

b_coal=solve(TPDM_hat[1:dim(Dat_coal_tr)[2]-1,1:dim(Dat_coal_tr)[2]-1])%*%TPDM_hat[1:dim(Dat_coal_tr)[2]-1,dim(Dat_coal_tr)[2]]
b_beer=solve(TPDM_hat[1:dim(Dat_coal_tr)[2]-1,1:dim(Dat_coal_tr)[2]-1])%*%TPDM_hat[1:dim(Dat_coal_tr)[2]-1,dim(Dat_coal_tr)[2]]
b_paper=solve(TPDM_hat[1:dim(Dat_coal_tr)[2]-1,1:dim(Dat_coal_tr)[2]-1])%*%TPDM_hat[1:dim(Dat_coal_tr)[2]-1,dim(Dat_coal_tr)[2]]

cbind.data.frame(b_coal,ColName_coal[1:29])
cbind.data.frame(b_beer,ColName_beer[1:29])
cbind.data.frame(b_paper,ColName_paper[1:29])

####  The best linear predictor
Xhat=Amul(t(b_coal),t(Dat_coal_te[,1:dim(Dat_coal_te)[2]-1]))
Xhat=as.vector(Xhat)
Xhat=Amul(t(b_beer),t(Dat_beer_te[,1:dim(Dat_coal_te)[2]-1]))
Xhat=as.vector(Xhat)
Xhat=Amul(t(b_paper),t(Dat_paper_te[,1:dim(Dat_coal_te)[2]-1]))
Xhat=as.vector(Xhat)
summary(Xhat)

####  Decompose the TPDM_Pred or the TPDM_Pred_hat
library(MASS)
source("C:/Linear_prediction/7_Ang_CPD.R")
TPDM_P_hat=TPDM_Out$TPDM_P_hat
AngCPD_Out=Ang_CPD(TPDM_P_hat = TPDM_P_hat,m = 80)

save(AngCPD_Out,file="C:/Linear_prediction/AngCPD_coal.Rdata")
save(AngCPD_Out,file="C:/Linear_prediction/AngCPD_beer.Rdata")
save(AngCPD_Out,file="C:/Linear_prediction/AngCPD_paper.Rdata")
load(file="C:/Linear_prediction/AngCPD_coal.Rdata")
load(file="C:/Linear_prediction/AngCPD_beer.Rdata")
load(file="C:/Linear_prediction/AngCPD_paper.Rdata")

####  Rescale transformed data
library(Hmisc)
####  Check the coverage rate
library(plotrix)
source("C:/Linear_prediction/8_JointRegion.R")
Out1=JointRegion(X_f = Dat_coal_te[,dim(Dat_coal_tr)[2]],Xhat = Xhat,Ang = AngCPD_Out$All_ang,Ang_mass = AngCPD_Out$All_mass,tol=0.001,plot = T,ray = FALSE,cone = T,xr1 = 0,xr2 = max(Xhat)*1.1,yr1 = 0,yr2 = max(Xhat)*1.1)
Out1=JointRegion(X_f = Dat_beer_te[,dim(Dat_coal_tr)[2]],Xhat = Xhat,Ang = AngCPD_Out$All_ang,Ang_mass = AngCPD_Out$All_mass,tol=0.001,plot = T,ray = FALSE,cone = T,xr1 = 0,xr2 = max(Xhat)*1.1,yr1 = 0,yr2 = max(Xhat)*1.1)
Out1=JointRegion(X_f = Dat_paper_te[,dim(Dat_coal_tr)[2]],Xhat = Xhat,Ang = AngCPD_Out$All_ang,Ang_mass = AngCPD_Out$All_mass,tol=0.001,plot = T,ray = FALSE,cone = T,xr1 = 0,xr2 = max(Xhat)*1.1,yr1 = 0,yr2 = max(Xhat)*1.1)
Out1$Coverage
#0.943
#0.883
#0.943

####  Plot of KDE for angular densities
library(VGAM)
library(ks)
source("C:/Linear_prediction/9_Ang_kde.R")
kde_cpd_out=Ang_kde(Ang = AngCPD_Out$All_ang, Ang_mass = AngCPD_Out$All_mass,bw = T,h = 0.3 ,plot_theta = T, plot_w = T,y_lim = 2)

####  Find the conditional interval
source("C:/Linear_prediction/10_CondDens.R")
a=2
CondDen_Out=CondDens(z = Out1$Mtx_P[a,],h_w2_CPD = kde_cpd_out$kde_trans_w, plot = T)
##  Inputs: cbind(X_predictant, Xhat), true kde, kde via CPD
##  Outputs: X_(p+1), conditional density, CDF / X_(p+1), conditional density via CPD, CDF via CPD

####  Find the prediction interval
source("C:/Linear_prediction/11_CondInterval.R")
Pred_inter=CondInterval(z2 = CondDen_Out$z2_CPD,X_f_single = Out1$Mtx_P[a,2],cumTraps = CondDen_Out$cumTraps,tol = 0.001)

####  Assess the coverage rate
source("C:/Linear_prediction/13_AssessCoverage.R")
AssessCoverage(Mtx_P = Out1$Mtx_P,h_w2_CPD = kde_cpd_out$kde_trans_w,Thres = 0.95,Thres_Xhat = T)
#Test_set[,1]
#W <- Test_set/Rad
#Keep <- Rad > quantile(Rad, .95)  
#Keep[Keep==TRUE]
#W_top <- W[Keep,]
#Mtx_P_top=Test_set[Keep,]

#0.979 (coverage rate (coal))
#0.963
#0.980

####  All conditional intervals for the top 5%
source("C:/Linear_prediction/14_CondIntervalPlot.R")
Cond_Out=CondIntervalPlot(Mtx_P = Out1$Mtx_P,Thres = 0.95,h_w2_CPD = kde_cpd_out$kde_trans_w)
Cond_mtx=Cond_Out$CondIntervalSave

####  Transform back to original scales
####  Coal
Fn=ECDF(negPf[11941:17911,28])
CDF=Fn(negPf[11941:17911,28])
CDF_c=cbind(negPf[11941:17911,28],CDF)
sl=sort.list(negPf[11941:17911,28])
CDF_sort=CDF_c[sl,]
CDF_sort=unique(CDF_sort)
#plot(x = CDF_sort[,1],y = CDF_sort[,2])
####  Beer
Fn=ECDF(negPf[11941:17911,29])
CDF=Fn(negPf[11941:17911,29])
CDF_c=cbind(negPf[11941:17911,29],CDF)
sl=sort.list(negPf[11941:17911,29])
CDF_sort=CDF_c[sl,]
CDF_sort=unique(CDF_sort)
####  Paper
Fn=ECDF(negPf[11941:17911,30])
CDF=Fn(negPf[11941:17911,30])
CDF_c=cbind(negPf[11941:17911,30],CDF)
sl=sort.list(negPf[11941:17911,30])
CDF_sort=CDF_c[sl,]
CDF_sort=unique(CDF_sort)

Xhat_orig=approxExtrap(x = CDF_sort[,2],y = CDF_sort[,1],xout = 1-(Mtx_P_top[,1]+shift)^(-2))
X_f_orig=approxExtrap(x = CDF_sort[,2],y = CDF_sort[,1],xout = 1-(Mtx_P_top[,2]+shift)^(-2))

Mtx_orig=cbind(Xhat_orig$y,X_f_orig$y)

##  Transform the scale of the conditional intervals back to the original scale
library(Hmisc)
#L_orig=approx(x = CDF_sort[,2],y = CDF_sort[,1],xout = 1-(Cond_intervals[,1]+shift)^(-2))
#U_orig=approx(x = CDF_sort[,2],y = CDF_sort[,1],xout = 1-(Cond_intervals[,2]+shift)^(-2))

L_orig=approxExtrap(x = CDF_sort[,2],y = CDF_sort[,1],xout = 1-(Cond_intervals[,1]+shift)^(-2))
U_orig=approxExtrap(x = CDF_sort[,2],y = CDF_sort[,1],xout = 1-(Cond_intervals[,2]+shift)^(-2))

L_o=(L_orig$y)
U_o=(U_orig$y)

####  All conditional intervals for the top 5%
dev.new()
par(mar=c(5.1,5.1,2,2))
#plot(Mtx_orig[,1],Mtx_orig[,2],xlim=c(0,max(Mtx_orig[,2])*1.5),ylim=c(0,max(Mtx_orig[,2])*1.5),
#     main=expression(paste("The scatter plot of"," ",X[p+1]," ","and"," ",hat(X)[p+1]," ","with conditional intervals")),xlab=expression(hat(X)[p+1]),ylab=expression(X[p+1]), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,pch=20,cex=1)
plot(Mtx_orig[,1],Mtx_orig[,2],xlim=c(0,max(Mtx_orig[,2])*1.5),ylim=c(0,max(Mtx_orig[,2])*1.5),
     main="",xlab="Predicted daily return",ylab="Observed daily return", cex.main=1.5, cex.lab=1.5, cex.axis=1.5,pch=20,cex=1.5)

for(a in 1:length(Mtx_orig[,1])){
  points(x = Mtx_orig[a,1],y=L_o[a],pch="-", col="blue",cex=2,lwd=2)
  points(x = Mtx_orig[a,1],y=U_o[a],pch="-", col="blue",cex=2,lwd=2)
}

####  Check the coverage rate with the original scale
sum(L_o <= Mtx_orig[,2] & Mtx_orig[,2] <= U_o)/length(Mtx_orig[,2])
#[1] 0.9186047    ## Same
