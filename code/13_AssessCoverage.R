AssessCoverage <- function(Mtx_P, h_w2_CPD, Thres=0.95, Thres_Xhat=T, Thres_rad=FALSE){
  
  ##  Choose the largest values given a large value of Xhat
  if(Thres_Xhat){
    Rad <- Mtx_P[,1]
    W <- Mtx_P/Rad
    Keep <- Rad > quantile(Rad, Thres)
    Keep[Keep==TRUE]
    W_top <- W[Keep,]
    Mtx_P_top=Mtx_P[Keep,]
  }
  ##  Choose the largest values given a norm is large
  if(Thres_rad){
    Rad <- sqrt(apply(Mtx_P^2, 1, sum)) 
    W <- Mtx_P/Rad
    Keep <- Rad > quantile(Rad, Thres)  
    Keep[Keep==TRUE]
    W_top <- W[Keep,]
    Mtx_P_top=Mtx_P[Keep,]
  }
  
  ##  Assess the coverage rate for conditonal density
  Counting=rep(0,dim(Mtx_P_top)[1])
  Pred_interval=matrix(0,nrow = dim(Mtx_P_top)[1],ncol = 2,byrow = T)
  for(i in 1:dim(Mtx_P_top)[1]){
    z2_ext=CondDens(z = Mtx_P_top[i,], h_w2_CPD = h_w2_CPD, plot = FALSE)  
    Count_ext=CondInterval(z2 = z2_ext$z2_CPD, X_f_single = Mtx_P_top[i,2], cumTraps = z2_ext$cumTraps, tol = 0.0001)
    Counting[i]=Count_ext$count
    Pred_interval[i,]=c(Count_ext$Con_L,Count_ext$Con_U)
    print(i)
    CoverageRate=sum(Counting)/dim(Mtx_P_top)[1]
  }
  
  return(list("CoverageRate"=CoverageRate,"Mtx_P_top"=Mtx_P_top))
}
