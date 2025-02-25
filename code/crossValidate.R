#' Cross validated bandwidth for angular density function
#' 
#' @description
#' Returns k-fold cross-validated bandwidth for kernel density estimate of 
#' the angular density.
#' 
#'
#' @param Dat Numeric data matrix
#' @param Ang Numeric vector of angular components
#' @param pMass Numeric vector of point masses
#' @param Thres Numeric of a high quantile for estimating TPDM
#' @param kfold Numeric of the k-fold
#' @param bandW Numeric of a bandwidth
#' @param Quan Numeric of a target high quantile
#'
#' @return Numeric of the average coverage rate
#' @export
#'
#' @examples 
#' ##  Optimal bandwidth using cross-validation
#' #targetQuan=0.97
#' #seqbw=seq(0,0.95,by=0.05)
#' #cv_result=sapply(seqbw,function(bw) crossValidate(Out=Out,kfold = 10,bandW = bw,Quan = targetQuan))
#' #cv_result
#' ##  Find the best bandwidth
#' #which.min(abs(cv_result-targetQuan))
crossValidate <- function(Dat, Ang, pMass, Thres, kfold=3, bandW, Quan){
  
  n=nrow(Dat);p=dim(Dat)[2]
  fold_size=floor(n/kfold)
  
  coverRate=c()
  
  for(i in 1:kfold){
    test_ind=((i-1)*fold_size+1):(i*fold_size)
    if(i==kfold){
      test_ind=((i-1)*fold_size+1):n
    }
    
    Train=Dat[-test_ind,]
    Test=Dat[test_ind,]
    
    Est=estimateParams(X = Train, Thres = Thres)
    b_vec=Est$bhat
    
    Xhat_test=Amul(t(b_vec),t(Test[,-p]))
    Xhat_test=as.vector(Xhat_test)
    jointOut=jointRegion(Xhat = Xhat_test, Xf = Test[,p],
                         Angular = Ang, Pmass = pMass, Quan = Quan,
                         Plot = F, axisLimit = 40, dataPoint = 1)
    
    if(bandW==0){
      kde_out=KDE_w(Ang = Ang,Pmass = pMass,Plot=F)
    }else{
      kde_out=KDE_w(Ang = Ang,Pmass = pMass,bw = T,h = bandW,Plot=F)
    }
    #coverOut=coverageRate(XY = jointOut$XY_top, kde_est = kde_out, Quan = Quan,Plot = F)
    
    XhatXp1 <- cbind(Xhat_test,Test[,p])
    Keep <- XhatXp1[,1] > quantile(XhatXp1[,1], Quan)
    XhatXp1_top <- XhatXp1[Keep,]
    coverOut=coverageRate(XY = XhatXp1_top, kde_est = kde_out, Quan = Quan,Plot = F)
    coverRate<-c(coverRate,coverOut$CoverageRate)
  }
  ##  Average coverage rate across folds
  mean(coverRate)
}


