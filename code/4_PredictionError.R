PredictionError <- function(Ncol=10,Nrow=4,n=20000,min=0,max=5,plot_scatt=FALSE,plot_TM=FALSE,plot_D=T){
  
  ##  Simulate a shifted Pareto dist
  U <- runif(n*Ncol)
  shift <- 0.9352074  # make the mean of InvT(Z) centered.
  Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol)
  
  ##  Generate a standardized matrix A from a uniform dist
  B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
  B_norm=sqrt(apply(B^2,1,sum))
  A <- B/B_norm
  
  ##  A known TPDM
  TPDM_X=A%*%t(A)
  
  ##  Tail ratio of prediction errors
  kstar=TPDM_X[Nrow,Nrow]-t(TPDM_X[1,2:Nrow])%*%solve(TPDM_X[1:(Nrow-1),1:(Nrow-1)])%*%TPDM_X[1,2:Nrow]
  
  ##  Simulate a random vector X by matrix multiplication
  X_t <- t(Amul(A, t(Z)))
  
  ##  Find the optimized vector b
  b=(solve(A[1:(Nrow-1),]%*%t(A[1:(Nrow-1),]),tol = 1e-30)%*%A[1:(Nrow-1),]%*%A[Nrow,]) # A[Nrow,] is a qx1 vector, not a 1xq.

  ##  The best transformed-linear predictor using TPDM
  Xhat=Amul(t(b),t(X_t[,1:(Nrow-1)]))
  Xhat=as.vector(Xhat)
  
  ##  Find the transformed-prediction error between Xhat and X_4
  Diff1=vSub(v1 = X_t[,Nrow],v2 = Xhat)
  Diff2=vSub(v1 = Xhat,v2 = X_t[,Nrow])
  D_TPDM=apply(X = cbind(Diff1,Diff2),1,max)

  ##  A scatterplot of "D" against Xhat with the 0.95 quantile for "D"
  if(plot_D){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(Xhat,D_TPDM,xlim=c(0,max(Xhat)),ylim=c(0,max(Xhat)),main="",xlab=expression(hat(X)[4]),ylab="D", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #plot(Xhat,D_TPDM,xlim=c(0,max(Xhat)),ylim=c(0,65),main="",xlab=expression(hat(X)[4]),ylab="D", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h=sqrt(2*kstar/0.05),lty=2,lwd=2)
    abline(h=sqrt(2*kstar/0.05),lty=2,lwd=2)
  }
  
  ##  Assess the coverage rate
  D_up=as.numeric(sqrt(kstar/0.025))
  Coverage_D=1-2*length(D_TPDM[D_TPDM>D_up])/length(D_TPDM)
  
  return(list("kstar"=kstar,"b"=b,"D_TPDM"=D_TPDM,"Coverage_D"=Coverage_D))
}


