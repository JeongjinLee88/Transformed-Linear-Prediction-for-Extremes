CondInterval <- function(z2, X_f_single, cumTraps, tol){
  
  ##  Sort the CDF
  Con_comb=cbind(z2,cumTraps)
  sl=sort.list(z2)
  Con_sort=Con_comb[sl,]
  
  ##  Interpolate the empirical CDF given a large value
  #z2[400]
  #z2_seq=seq(min(z2),z2[400],length.out = 10000)
  z2_seq=c(seq(0,10,by=0.00001),seq(10.001,300,length.out=100000))
  Int_Con=approx(x = Con_sort[,1], y = Con_sort[,2], xout = z2_seq, method = "linear")
  Loc_Con_L=which(x = Int_Con$y-0.025 > 0 & Int_Con$y-0.025 < tol)
  Loc_Con_U=which(x = Int_Con$y-0.975 > 0 & Int_Con$y-0.975 < tol)
  
  ##  Find bounds
  Con_L=max(Int_Con$x[Loc_Con_L])
  Con_U=min(Int_Con$x[Loc_Con_U])
  
  X_future=X_f_single
  if(Con_L <= X_future & X_future <= Con_U){
    count=1
  } else{
    count=0
  }
  
  return(list("Con_L"=Con_L, "Con_U"=Con_U, "count"=count))
}

