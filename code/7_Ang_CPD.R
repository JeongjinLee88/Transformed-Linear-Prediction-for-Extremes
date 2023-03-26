Ang_CPD <- function(TPDM_P_hat=TPDM_P_hat, Add_c=8, m=80){
  
  ##  TPDM_Pred: CP and rank=2
  ##  An initial factorization TPDM_Pred=BB' (Choleski)
  B=t(chol(x = TPDM_P_hat)) 
  #t(chol(x = TPDM_P_hat))%*%chol(x = TPDM_P_hat)
  
  ##  Append columns to B
  B_app=cbind(B[,2],matrix(rep(B[,1],Add_c)/sqrt(Add_c),nrow=2))
  #B_app%*%t(B_app)

  ##  CPD(Completely Positive Decomposition)
  BQ_save=array(0,dim=c(2,1+Add_c,m))
  for(j in 1:m){
    ##  Generate a rxr random matrix M
    M=matrix(rnorm((1+Add_c)*(1+Add_c),mean = 0,sd = 1),nrow=(1+Add_c))
    #M=matrix(runif((1+Add_c)*(1+Add_c)),nrow=(1+Add_c))
    ##  Find Q_0 by projecting M onto Q_r
    Out1=eigen(M%*%t(M))
    U=Out1$vectors
    Out2=eigen(t(M)%*%M)
    V=Out2$vectors
    ##  Initial orthogonal matrix
    Q_0=U%*%t(V)
    ##  Iterate
    BQ=B_app%*%Q_0
    for(i in 1:5000){
      BQ[BQ<=0]<-0
      D=BQ
      B_plus=ginv(X = D)
      Phat=B_plus%*%D+(diag((1+Add_c))-B_plus%*%B_app)%*%Q_0
      Out11=eigen(Phat%*%t(Phat))
      Out22=eigen(t(Phat)%*%Phat)
      Q_0=Out11$vectors%*%t(Out22$vectors)
      BQ=B_app%*%Q_0
      if(all(BQ>=-10^(-15))){
        break
      }
      print(i)
    }
    BQ_save[,,j]=BQ
  }

  for(k in 1:m){
    if(all(BQ_save[,,k]>=0)){
      BQ_save[,,k]=BQ_save[,,k]  
    }
    else{
      BQ_save[,,k]=0
    }
  }
  
  Nonzero=BQ_save[BQ_save!=array(0,dim=c(2,(1+Add_c),m))]
  Nonzero=data.frame(Nonzero)
  m_star=length(unlist(Nonzero))/(2*(1+Add_c))
  BQ_array=array(unlist(Nonzero),dim=c(2,(1+Add_c),m_star))
  
  ##  Find estimated angular components under repeated decompositions
  All_ang=matrix(rep(0,(1+Add_c)*m_star),nrow = m_star)
  All_mass=matrix(rep(0,(1+Add_c)*m_star),nrow = m_star)
  for(j in 1:m_star){
    for(i in 1:(1+Add_c)){
      All_ang[j,i]=atan(x = BQ_array[2,i,j]/BQ_array[1,i,j])
      All_mass[j,i]=sum(BQ_array[,i,j]^2)/m_star
    }
  }
  
  All_ang=as.vector(All_ang)
  All_mass=as.vector(All_mass)
  
  return(list("All_ang"=All_ang,"All_mass"=All_mass))  
}





