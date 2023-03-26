TPDM_Ang <- function(Ncol=400, Nrow=7, n=60000, min=0, max=5){
  
  ##  Simulate a shifted Pareto dist
  U <- runif(n*Ncol)
  shift <- 0.9352074  # make the mean of InvT(Z) be 0.
  Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol)
  
  ##  Generate a matrix A
  B <- matrix(runif(Nrow*Ncol,min = min,max = max), nrow = Nrow, ncol = Ncol)
  B_norm=sqrt(apply(B^2,1,sum))
  A <- B/B_norm  # Normalized
  
  ##  TPDM
  TPDM_X=A%*%t(A)  
  
  ##  Generate a vector X by matrix multiplication
  X_t <- t(Amul(A, t(Z)))
  
  ##  Generate a vector b
  b=solve(TPDM_X[1:(Nrow-1),1:(Nrow-1)])%*%TPDM_X[1:(Nrow-1),Nrow]
  
  ##  The best transformed-linear predictor
  Xhat=Amul(t(b),t(X_t[,1:(Nrow-1)]))
  Xhat=as.vector(Xhat)
  
  ##  The predictor TPDM of Xhat and X
  A_Pred=rbind(t(b)%*%A[1:(Nrow-1),],A[Nrow,])  #2xq
  TPDM_Pred=A_Pred%*%t(A_Pred)
  
  ##  Scaled A_Pred to have a total mass 1
  A_S=A_Pred/sqrt(sum(diag(TPDM_Pred)))
  TPDM_Prob=A_S%*%t(A_S)
  
  ##  The squared scale of D
  TR=TPDM_X[Nrow,Nrow]-t(TPDM_X[1,2:Nrow])%*%solve(TPDM_X[1:(Nrow-1),1:(Nrow-1)])%*%TPDM_X[1,2:Nrow]
  TR
  
  ##  Known angular components and masses
  Ang_T=atan(A_Pred[2,]/A_Pred[1,])
  Ang_mass=apply(A_Pred^2,2,sum)
  
  
  return(list("X_t"=X_t, "Xhat"=Xhat, "b"=b, "A_Pred"=A_Pred, "TPDM_Pred"=TPDM_Pred, "TR"=TR, "Ang_T"=Ang_T, "Ang_mass"=Ang_mass, "A"=A))  
}



