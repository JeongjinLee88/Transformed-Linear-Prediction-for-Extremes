#' Pairwise estimation for extremal dependence measure
#' @description
#' Perform pairwise estimation for extremal dependence measure using high threshold 
#' of \eqn{\|\widehat{X}_{p+1},X_{p+1}\|_2}.
#' 
#' @param X \eqn{n\times p} numeric data matrix with Pareto variables
#' @param Thres Numeric; a high threshold quantile
#'
#' @return a list of numeric results:
#' 1. TPDM estimate 2. prediction TPDM estimate 3. the estimated vector \eqn{\boldsymbol{b}}
#' @export
#'
#' @examples
estimateParams <- function(X,Thres=0.95){
  
  d=dim(X)[2]  
  
  ## Pairwise estimation
  TPDM_hat=matrix(0,d,d)
  for(i in 1:d){
    for(j in 1:i){
      R_pair=sqrt(apply(X[,c(i,j)]^2,1,sum))
      W_pair=X[,c(i,j)]/R_pair
      Keep=R_pair > quantile(R_pair, Thres)
      W_keep=W_pair[Keep,]
      TPDM_hat[i,j]=t(W_keep[,1])%*%W_keep[,2] / dim(W_keep)[1] * 2
    }
  }
  ##  Estimated TPDM
  TPDM_hat=TPDM_hat+t(TPDM_hat)-diag(diag(TPDM_hat))
  
  Sig11=TPDM_hat[1:(d-1),1:(d-1)]
  Sig12=TPDM_hat[1:(d-1),d]
  Sig21=t(Sig12)
  Sig22=TPDM_hat[d,d]
  
  ##  Prediction inner product matrix
  Sig_P_11=Sig21%*%solve(Sig11)%*%Sig12
  TPDM_Phat=matrix(c(Sig_P_11,Sig_P_11,Sig_P_11,Sig22),nrow=2)
  ##  b hat
  bhat=solve(Sig11)%*%Sig12
  
  return(list("TPDM_hat"=TPDM_hat,"TPDM_Phat"=TPDM_Phat,"bhat"=bhat))
}
