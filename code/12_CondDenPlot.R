CondDenPlot <- function(a, X_f, Xhat, Ang, Ang_mass, h_w2_CPD = kde_cpd_out$kde_trans_w, tol=0.0001, ConDen_CPD=T, Cond_Int_CPD=T, xr1=25.84, xr2=25.86, yr1=0, yr2=50){
  
  ##  Plot a joint polar region
  par(mar=c(5.1,5.1,2,2))
  Out1=JointRegion(X_f = X_f, Xhat = Xhat, Ang = Ang, Ang_mass = Ang_mass, tol=tol, plot = T, ray = FALSE, cone=FALSE, xr1=xr1, xr2=xr2, yr1=yr1, yr2=yr2)
  
  for(a in a){
    
    ##  The approximate conditional density given a large observed value
    CondDen_Out=CondDens(z = Out1$Mtx_P[a,], h_w2_CPD = h_w2_CPD, plot = FALSE)
    
    if(ConDen_CPD){
      lines(Out1$Mtx_P[a,1]+CondDen_Out$densVal, CondDen_Out$z2_CPD, lwd=2)    
    }
    
    ##  Indicate a large observed value
    abline(v = Out1$Mtx_P[a,1])
    text(x = Out1$Mtx_P[a,1]+0.3,y = -0.2,labels = round(Out1$Mtx_P[a,1],2))
    points(x = Out1$Mtx_P[a,1], y=Out1$Mtx_P[a,2], pch="*", col="blue", cex=3, lwd=2)
    
    if(Cond_Int_CPD){
      ##  Conditional prediction intervals
      Pred_inter_CPD=CondInterval(z2 = CondDen_Out$z2_CPD, X_f_single = Out1$Mtx_P[a,2], cumTraps = CondDen_Out$cumTraps, tol = 0.0001)
      points(x = Out1$Mtx_P[a,1],y=Pred_inter_CPD$Con_L,pch="-", col="red",cex=3,lwd=2)
      points(x = Out1$Mtx_P[a,1],y=Pred_inter_CPD$Con_U,pch="-", col="red",cex=3,lwd=2)
      print(c(Pred_inter_CPD$Con_L,Pred_inter_CPD$Con_U))
      print(Out1$Mtx_P[a,])
    }
  }
}

