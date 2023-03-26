CondIntervalPlot <- function(Mtx_P, Thres, h_w2_CPD){
  
  ##  Choose the largest values given Xhat is large
  Rad <- Mtx_P[,1]
  W <- Mtx_P/Rad
  Keep <- Rad > quantile(Rad, Thres)
  Keep[Keep==TRUE]
  W_top <- W[Keep,]
  Mtx_top=Mtx_P[Keep,]
  
  ##  A scatterplot of X_p against Xhat given Xhar is large
  dev.new()
  par(mar=c(5.1,5.1,2,2))
  plot(Mtx_top[,1], Mtx_top[,2], xlim=c(0,20), ylim=c(0,35)
       , main="", xlab="", ylab=expression(X[p+1]), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3)
  title(xlab=expression(hat(X)[p+1]),line=3.5, cex.lab=1.5)
  Qt=quantile(Rad, Thres)
  abline(v=Qt,lty=2,lwd=2)
  
  ##  Save conditional intervals
  CondIntervalSave=matrix(rep(NA,length(Mtx_top[,1])*4),nrow = length(Mtx_top[,1]),ncol = 4)
  for(a in 1:length(Mtx_top[,1])){
    ##  Find condtional densities
    CondDen_Out=CondDens(z = Mtx_top[a,], h_w2_CPD = h_w2_CPD, plot = FALSE)
    ##  Find conditional intervals
    CondInterval_Out=CondInterval(z2 = CondDen_Out$z2, X_f_single = Mtx_top[a,2], cumTraps = CondDen_Out$cumTraps, tol = 0.0001)
    ##  Draw conditoinal intervals
    points(x = Mtx_top[a,1],y=CondInterval_Out$Con_L,pch="-", col="blue",cex=2,lwd=2)
    points(x = Mtx_top[a,1],y=CondInterval_Out$Con_U,pch="-", col="blue",cex=2,lwd=2)
    print(c(CondInterval_Out$Con_L,CondInterval_Out$Con_U,Mtx_top[a,2]))
    CondIntervalSave[a,]=c(CondInterval_Out$Con_L,CondInterval_Out$Con_U,Mtx_top[a,])
  }
  return(list("CondIntervalSave"=CondIntervalSave))
}
