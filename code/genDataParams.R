#' Simulate from regularly varying random vector using the true matrix A
#' 
#' @description
#' `genDataParams()` creates a \eqn{n \times p} data matrix \eqn{\boldsymbol{X}}
#' constructed from a matrix \eqn{A} applied to a random vector \eqn{\boldsymbol{Z}} of independent unit-scale
#' Pareto random variables with \eqn{\alpha=2}.
#' The function computes relevant true quantities such as the optimal vector b, true TPDM, prediction TPDM, 
#' prediction error, true angular components, and true angular masses.
#' 
#' @param A \eqn{p \times q} numeric matrix; set q > p. 
#' @param Z \eqn{n \times q} numeric matrix; a unit-scale Pareto distribution with \eqn{\alpha=2}. 
#'
#' @return a list of generated data with relevant true quantities
#' 
#' @export
#'
#' @examples 
#' U <- runif(n*Ncol)
#' shift <- 0.9352074  # make the mean of InvT(Z) centered.
#' Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol)
#' B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
#' B_norm=sqrt(apply(B^2,1,sum))
#' A <- B/B_norm
#' genDataParams(A,Z)
genDataParams <- function(A, Z){
  
  p=dim(A)[1]
  ##  Specified TPDM
  TPDM_X=A%*%t(A)
  
  ##  Generate a vector X by matrix multiplication
  Xp <- t(Amul(A, t(Z)))
  
  ##  Optimal vector b = Sig11^(-1) * Sig12
  b=solve(TPDM_X[1:(p-1),1:(p-1)])%*%TPDM_X[1:(p-1),p] 
  
  ##  The best transformed-linear predictor
  Xhat=Amul(t(b),t(Xp[,1:(p-1)]))
  Xhat=as.vector(Xhat)
  
  ##  The predictor TPDM of Xhat and X
  A_Pred=rbind(t(b)%*%A[1:(p-1),],A[p,])  #2xq
  TPDM_Pred=A_Pred%*%t(A_Pred)
  
  ##  The squared scale of D
  PredE=TPDM_X[p,p]-TPDM_X[p,1:(p-1)]%*%solve(TPDM_X[1:(p-1),1:(p-1)])%*%TPDM_X[p,1:(p-1)]
  PredE
  
  ##  Known angular components and masses
  Ang_T=atan(A_Pred[2,]/A_Pred[1,])
  Ang_mass=apply(A_Pred^2,2,sum)
  
  return(list("Xp"=Xp, "TPDM_X"=TPDM_X, "Xhat"=Xhat, "b"=b, "A_Pred"=A_Pred, "TPDM_Pred"=TPDM_Pred, "PredE"=PredE, "Ang_T"=Ang_T, "Ang_mass"=Ang_mass))  
}



