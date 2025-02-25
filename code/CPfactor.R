#' Completely positive factorization for a matrix
#' 
#' @description
#' `CPfactor()` conducts a completely positive factorization to the \eqn{2\times 2}
#' prediction (estimated) TPDM.
#' This function can be applied to a higher dimensional matrix in principle.
#' 
#' 
#' @param Mtx Numeric square matrix
#' @param q_star Numeric; the column dimension of a matrix A
#' @param ite_pereach Numeric value; the number of iterations until the algorithm converges
#' @param ite_cp Numeric value; the number of CP-factorization
#'
#' @return a list of two numeric vectors; `angular` is a vector of angular components
#' and `pmass` is a vector of point masses derived from the matrix \eqn{A} 
#' @export
#'
#' @examples
#' TPDM_Phat=matrix(c(rep(0.7540613,3),1),2)
#' CPfactor(TPDM_Phat, q_star=9, ite_pereach=5000, ite_cp=5)
CPfactor <- function(Mtx, q_star=9, ite_pereach=5000, ite_cp=80){
  
  dcol=dim(Mtx)[2]
  if(q_star < dcol){
    stop("q_star must be greater than  or equal to the column dimension of the given matrix")
  }
  add_col=q_star-(dcol-1)
  
  ##  TPDM_Pred (Rank=2)
  if(any(eigen(x = Mtx)$values < 0)){
    Mtx=nearPD(Mtx)
  }
  
  ##  Cholesky factorization for TPDM_Pred=BB' (Initial factorization)
  Btilde=t(chol(Mtx)) 
  #t(chol(x = Mtx))%*%chol(x = Mtx) = TPDM_Pred
  
  ##  Add extra columns to B
  AugmentedMatrix=matrix(rep(Btilde[,1],add_col)/sqrt(add_col),nrow=dim(Btilde)[1])
  B=cbind(Btilde[,2],AugmentedMatrix)
  #B%*%t(B) = TPDM_Pred

  ##  CP-factorization
  BQ_save=array(0,dim=c(dcol,q_star,ite_cp))
  for(j in 1:ite_cp){
    converged <- FALSE
    while(!converged){
      ##  Generate a qxq random matrix M
      M=matrix(rnorm((q_star)*(q_star)),nrow=(q_star))
      ##  Find Q_0 by projecting M onto O_r
      U=eigen(M%*%t(M))$vectors
      V=eigen(t(M)%*%M)$vectors
      ##  Initial orthogonal matrix
      Q_0=U%*%t(V)
      
      ##  Alternating projections
      BQ=B%*%Q_0
      for(i in 1:ite_pereach){
        BQ[BQ<=0]<-0
        D=BQ
        PseInv=ginv(X = D)
        Phat=PseInv%*%D+(diag((q_star))-PseInv%*%B)%*%Q_0
        LeftSingular=eigen(Phat%*%t(Phat))$vectors
        RightSingular=eigen(t(Phat)%*%Phat)$vectors
        Q_0=LeftSingular%*%t(RightSingular)
        BQ=B%*%Q_0
        BQ=Re(BQ)
        if(all(BQ >= -10^(-15))){
          converged <- TRUE
          break
        }
        #print(i)
      }
    }
    BQ_save[,,j]=BQ
  }

  array_reorder=aperm(BQ_save,c(1,3,2)) # reorder dimensions to (dcol x ite_cp x q_star)
  BQ_comb=matrix(array_reorder, nrow=dcol)
  
  ##  save angular points and point masses
  pmass=apply(BQ_comb^2,2,sum)/ite_cp
  angular=BQ_comb[1,]/sqrt(apply(BQ_comb^2,2,sum))
  
  return(list("angular"=angular,"pmass"=pmass))  
}





