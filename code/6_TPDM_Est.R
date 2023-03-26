TPDM_Est <- function(Nrow=7, Thres=0.99, X_t=X_t){
  
  ## Pairwise estimation
  TPDM_hat=matrix(0,nrow = Nrow,ncol = Nrow)
  for(i in 1:Nrow){
    for(j in 1:i){
      R_pair=sqrt(apply(X_t[,c(i,j)]^2,1,sum))
      W_pair=X_t[,c(i,j)]/R_pair
      Keep=R_pair > quantile(R_pair, Thres)
      W_keep=W_pair[Keep,]
      TPDM_hat[i,j]=t(W_keep[,1])%*%W_keep[,2] / dim(W_keep)[1] * 2
    }
  }
  ##  Estimated TPDM
  TPDM_hat=TPDM_hat+t(TPDM_hat)-diag(diag(TPDM_hat))
  
  Sig11=TPDM_hat[1:(Nrow-1),1:(Nrow-1)]
  Sig12=TPDM_hat[1:(Nrow-1),Nrow]
  Sig21=t(Sig12)
  Sig22=TPDM_hat[Nrow,Nrow]
  
  Sig_P_11=Sig21%*%solve(Sig11)%*%Sig12
  ##  Prediction inner product matrix
  TPDM_P_hat=matrix(c(Sig_P_11,Sig_P_11,Sig_P_11,Sig22),nrow=2)
  TPDM_P_hat
  bhat=solve(Sig11)%*%Sig12
  
  return(list("TPDM_hat"=TPDM_hat,"TPDM_P_hat"=TPDM_P_hat,"bhat"=bhat))
}
