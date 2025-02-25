####  Load portfolios data
setwd("/home/leej40/Documents/extlinear/data")
Pf <- read.csv("30_Industry_Portfolios_Daily.CSV")
length(Pf[,1]) #24854
Pf[Pf==-99.99]<-NA
Pf=na.omit(Pf)

####  Select data for 1950-2020
Pf <- Pf[Pf[,1] > 19500000,]
dates <- Pf[,1]
Pf <- Pf[, -1]
d <- dim(Pf)[2]
n <- dim(Pf)[1]  #17911

####  Reorder variables (last columns: coal, beer, paper)
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

save(negPf_train,file="negPf_train.Rdata")
save(negPf_test,file="negPf_test.Rdata")
save(N_train,file="N_train.Rdata")
save(N_test,file="N_test.Rdata")

load(file="negPf_train.Rdata")
load(file="negPf_test.Rdata")

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
#0.7926421
##  Beer
Comb_n=cbind(Xhat_n_beer,N_beer_te[,30])
Q_n=quantile(x = Xhat_n_beer,probs = 0.95)
Top_n=Comb_n[Comb_n[,1]>Q_n,]

UB=Top_n[,1]+1.96*sqrt(MSPE_beer)
LB=Top_n[,1]-1.96*sqrt(MSPE_beer)

Count_n=Top_n[,2]>LB & Top_n[,2]<UB
sum(Count_n)/length(Count_n)
#0.8394649
##  Paper
Comb_n=cbind(Xhat_n_paper,N_paper_te[,30])
Q_n=quantile(x = Xhat_n_paper,probs = 0.95)
Top_n=Comb_n[Comb_n[,1]>Q_n,]

UB=Top_n[,1]+1.96*sqrt(MSPE_paper)
LB=Top_n[,1]-1.96*sqrt(MSPE_paper)

Count_n=Top_n[,2]>LB & Top_n[,2]<UB
sum(Count_n)/length(Count_n)
#0.9230769
####  Estimate the TPDM
source("estimateParams.R")
Thres_u=0.95
Est_coal=estimateParams(X = Dat_coal_tr, Thres = Thres_u) # Pareto variables
Est_beer=estimateParams(X = Dat_beer_tr, Thres = Thres_u) # Pareto variables
Est_paper=estimateParams(X = Dat_paper_tr, Thres = Thres_u) # Pareto variables

####  Find the optimized Vector b
eigen(Est_coal$TPDM_hat)$values
eigen(Est_beer$TPDM_hat)$values
eigen(Est_paper$TPDM_hat)$values

b_coal=Est_coal$bhat
b_beer=Est_beer$bhat
b_paper=Est_paper$bhat

#cbind.data.frame(b_coal,ColName_coal[1:29])
#cbind.data.frame(b_beer,ColName_beer[1:29])
#cbind.data.frame(b_paper,ColName_paper[1:29])

####  The best linear predictor
Xhat_coal=Amul(t(b_coal),t(Dat_coal_te[,1:dim(Dat_coal_te)[2]-1]))
Xhat_coal=as.vector(Xhat_coal)
Xhat_beer=Amul(t(b_beer),t(Dat_beer_te[,1:dim(Dat_coal_te)[2]-1]))
Xhat_beer=as.vector(Xhat_beer)
Xhat_paper=Amul(t(b_paper),t(Dat_paper_te[,1:dim(Dat_coal_te)[2]-1]))
Xhat_paper=as.vector(Xhat_paper)


####  Decompose the TPDM_Pred or the TPDM_Pred_hat
library(MASS)
source("CPfactor.R")
CP_coal=CPfactor(Mtx = Est_coal$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 50)
CP_beer=CPfactor(Mtx = Est_beer$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 50)
CP_paper=CPfactor(Mtx = Est_paper$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 50)
save(CP_coal,file="CP_coal.RData")
save(CP_beer,file="CP_beer.RData")
save(CP_paper,file="CP_paper.RData")
load(file="CP_coal.Rdata")
load(file="CP_beer.Rdata")
load(file="CP_paper.Rdata")

####  Check the coverage rate
####  Create the 95% joint polar region from a CP-factorization
library(plotrix)  # 'draw.circle'
source("jointRegion.R")
##  Find the 95% joint polar region
jointOut_coal=jointRegion(Xhat = Xhat_coal, Xf = Dat_coal_te[,30],
                     Angular = CP_coal$angular, Pmass = CP_coal$pmass, Quan = 0.95,
                     Plot = T, axisLimit = 40, dataPoint = 1)
jointOut_coal$coverage

jointOut_beer=jointRegion(Xhat = Xhat_beer, Xf = Dat_beer_te[,30],
                     Angular = CP_beer$angular, Pmass = CP_beer$pmass, Quan = 0.95,
                     Plot = T, axisLimit = 40, dataPoint = 1)
jointOut_beer$coverage

jointOut_paper=jointRegion(Xhat = Xhat_paper, Xf = Dat_paper_te[,30],
                     Angular = CP_paper$angular, Pmass = CP_paper$pmass, Quan = 0.95,
                     Plot = T, axisLimit = 40, dataPoint = 1)
jointOut_paper$coverage

##  Find the bandwidth via cross-validation
seqbw=seq(0.1,0.7,by=0.05) # a seq of bandwidth
cv_bw=rep(NA,length(seqQuan))
## coverage rates for the target quantile
cv_result=sapply(seqbw,function(bw) crossValidate(Dat = Dat_coal_te[,-30],Ang =CP_coal$angular,
                                                  pMass = CP_coal$pmass,Thres = Thres_u,
                                                  kfold = 3,
                                                  bandW = bw,Quan = 0.95))
cv_coal=seqbw[which.min(abs(cv_result-0.95))]

cv_result=sapply(seqbw,function(bw) crossValidate(Dat = Dat_beer_te[,-30],Ang =CP_beer$angular,
                                                  pMass = CP_beer$pmass,Thres = Thres_u,
                                                  kfold = 3,
                                                  bandW = bw,Quan = 0.95))
cv_beer=seqbw[which.min(abs(cv_result-0.95))]

cv_paper=sapply(seqbw,function(bw) crossValidate(Dat = Dat_paper_te[,-30],Ang =CP_paper$angular,
                                                  pMass = CP_paper$pmass,Thres = Thres_u,
                                                  kfold = 3,
                                                  bandW = bw,Quan = 0.95))
cv_paper=seqbw[which.min(abs(cv_result-0.95))]

####  Plot of KDE for angular densities
library(VGAM)
library(ks)
source("KDE_w.R")
kde_out_coal=KDE_w(Ang = CP_coal$angular,Pmass = CP_coal$pmass,bw = T,h = cv_bw, Plot=T)
kde_out_beer=KDE_w(Ang = CP_beer$angular,Pmass = CP_beer$pmass,bw = T,h = cv_bw, Plot=T)
kde_out_paper=KDE_w(Ang = CP_paper$angular,Pmass = CP_paper$pmass,bw = T,h = cv_bw, Plot=T)

kde_out_coal=KDE_w(Ang = CP_coal$angular,Pmass = CP_coal$pmass, Plot=T)
kde_out_beer=KDE_w(Ang = CP_beer$angular,Pmass = CP_beer$pmass, Plot=T)
kde_out_paper=KDE_w(Ang = CP_paper$angular,Pmass = CP_paper$pmass,Plot=T)

####  Assess the coverage rate
####  Plot conditional intervals with lines
source("coverageRate.R")
target_rate=0.95
XhatXp1 <- cbind(Xhat_coal,Dat_coal_te[,30])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out_coal, Quan = target_rate,Plot = T)
coverOut$CoverageRate # 0.9531773

XhatXp1 <- cbind(Xhat_beer,Dat_beer_te[,30])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out_beer, Quan = target_rate,Plot = T)
coverOut$CoverageRate # 0.9531773

XhatXp1 <- cbind(Xhat_paper,Dat_paper_te[,30])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out_paper, Quan = target_rate,Plot = T)
coverOut$CoverageRate # 0.9331104

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
