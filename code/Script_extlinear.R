##  Set the proper directory
#setwd("")
setwd("/home/leej40/Documents/extlinear/code")
##  Define transformed operations
source("TransformedOperations.R")

##  Calculate the tail ratio of prediction error
source("PredictionError.R")
set.seed(1205)
PredictionError(Ncol = 10, Nrow = 4, n = 20000, min = 0, max = 5, plot_D=T)

##  Simulate a shifted Pareto dist
set.seed(1234)
Nrow=7; Ncol=400; n=60000; min=0; max=5
U <- runif(n*Ncol)
shift <- 0.9352074  # make the mean of InvT(Z) centered.
Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol) # Necessary
##  Generate a p x q matrix A from a uniform dist
B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
B_norm=sqrt(apply(B^2,1,sum))
A <- B/B_norm

##  Generate X = A o Z
source("genDataParams.R")
Out=genDataParams(A = A, Z = Z)
n_train=ceiling(n*(2/3))
Train=Out$Xp[1:n_train,]
Test=Out$Xp[-(1:n_train),]

##  Plot of the approximate true angular measure
#pmass_Apred=apply(Out$A_Pred^2,2,sum)
#angular_Apred=Out$A_Pred[1,]/sqrt(apply(Out$A_Pred^2,2,sum))
#dev.new()
#plot(angular_Apred,pmass_Apred)

##  Estimate the TPDM, TPDM_pred, and b.
source("estimateParams.R")
Thres_u=0.75 # tried 0.9 and 0.95
Est=estimateParams(X = Train, Thres = Thres_u)
Est$TPDM_hat
Out$TPDM_X # true TPDM

Est$TPDM_Phat
Out$TPDM_Pred # true prediction TPDM

##  CP-factorization for a 2x2 TPDM_Phat
library(Matrix) #nearPD
library(MASS) #ginv
source("CPfactor.R")
CPout=CPfactor(Mtx = Est$TPDM_Phat,q_star = 10,ite_pereach = 5000,ite_cp = 100)
#CPout=CPfactor(Mtx = Out$TPDM_Pred,q_star = 10,ite_pereach = 5000,ite_cp = 100)
plot(CPout$angular,CPout$pmass)

##  Save output from TPDM_Phat
#save(CPout,file="CPout.RData")
load(file="CPout.RData")

##  Create the 95% joint polar region from a CP-factorization
library(plotrix)  # 'draw.circle'
source("jointRegion.R")
##  Calculate Xhat in the test set
Xhat_test=Amul(t(Est$bhat),t(Test[,-Nrow]))
Xhat_test=as.vector(Xhat_test)
##  Find the 95% joint polar region in Figure 2 (left)
target_rate=0.94
jointOut=jointRegion(Xhat = Xhat_test, Xf = Test[,Nrow],
                     Angular = CPout$angular, Pmass = CPout$pmass, Quan = target_rate,
                     Plot = T, axisLimit = 80, dataPoint = 471)
jointOut$coverage
# 471th obs corresponds to (Xhat,Xp+1)=(27.06545, 40.945)

##  Find the bandwidth via cross-validation for the target rate of 0.95
source("crossValidate.R")
target_rate=0.94
#Thres_u=0.75
seqbw=seq(0.1,0.7,by=0.05) # a seq of bandwidth
cv_result=sapply(seqbw,function(bw) crossValidate(Dat = Out$Xp,Ang =CPout$angular,
                                                  pMass = CPout$pmass,Thres = Thres_u,
                                                  kfold = 5,
                                                  bandW = bw,Quan = target_rate))
##  Find cv-bandwidths
cv_bw=seqbw[which.min(abs(cv_result-target_rate))]
cv_bw
# 0.4 for target rate of 0.90
# 0.35 for target rate of 0.95
# 0.6 for target rate of 0.98 

##  Kernel density estimation to approximate an angular density 'h'
##  Angular components and masses obtained from a CP-factor
library(VGAM) # probit
library(ks) # kde
source("KDE_w.R")
kde_out=KDE_w(Ang = CPout$angular, Pmass = CPout$pmass, bw = T, h = cv_bw, Plot=T)
kde_out

##  Plot an approximate conditional density with the 95% conditional interval
##  to reproduce Figure 2 (middle)
source("condDensity.R")
conden_out=condDensity(xhatPoint=Xhat_test[471],xp1Point=Test[471,Nrow],
                       kde_h = kde_out, Quan = target_rate,
                       xlim_upper = Xhat_test[471]+0.06, ylim_upper = 90, Plot_h = FALSE)

##  Assess the coverage rate
##  Plot conditional intervals using the cross-validated bandwidth
##  To reproduce Figure 2 (right)
source("coverageRate.R")
target_rate=0.94
XhatXp1 <- cbind(Xhat_test,Test[,Nrow])
Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], target_rate)
XhatXp1_top <- XhatXp1[Keep,]
coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = target_rate, Plot = T)
coverOut$CoverageRate
# 0.89
# 0.956
# 0.98
