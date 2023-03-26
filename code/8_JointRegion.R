JointRegion <- function(X_f, Xhat, Ang, Ang_mass, tol, plot=T, ray=T, cone=T,xr1=0,xr2=10,yr1=0,yr2=10,a=15171,Plot_Point=FALSE){
  
  ##  Empirical angular measure
  Ang_c=cbind(Ang,Ang_mass)
  sl=sort.list(Ang)
  Ang_sort=Ang_c[sl,]
  Cum_s=cumsum(Ang_sort[,2])
  Ang_ecdf=cbind(Ang_sort,Cum_s)
  
  ##  Total mass
  Total_mass=sum(Ang_mass)
  
  ##  Find lower and upper bounds at 0.025 and 0.975 empirical quantiles
  Ang_seq=seq(0,pi/2,length.out = 100000)
  Int_ang=approx(x=Ang_ecdf[,1], y = Ang_ecdf[,3], xout=Ang_seq, method = "linear")
  Loc_L=which(x = abs(Int_ang$y-Total_mass*0.025) < tol)
  Loc_U=which(x = abs(Int_ang$y-Total_mass*0.975) < tol)
  
  Ang_L=max(Int_ang$x[Loc_L])
  Ang_U=min(Int_ang$x[Loc_U])
  
  ##  Transforming to polar coordinates
  Mtx_P <- cbind(Xhat,X_f)
  Rad <- sqrt(apply(Mtx_P^2, 1, sum))
  W <- Mtx_P/Rad
  Keep <- Rad > quantile(Rad, .95)
  Keep[Keep==TRUE]
  W_top <- W[Keep,]
  Mtx_top <- Mtx_P[Keep,]
  
  ##  Assess a coverage rate in terms of the slopes of angular components
  Slp=W_top[,2]/W_top[,1]
  L_rate=length(Slp[Slp<tan(Ang_L)])/length(Slp)
  U_rate=length(Slp[Slp>tan(Ang_U)])/length(Slp)
  Coverage=1-(L_rate+U_rate)
  
  ##  A scatterplot of Xhat against X_f with a joint polar region
  X_comb=data.frame(Mtx_P,Rad)
  X_comb$color=as.character(cut(X_comb[,3],breaks = c(0,as.numeric(quantile(X_comb[,3],0.95)),Inf),labels = c("grey","black"),right = FALSE))
  R_95=quantile(Rad,probs = 0.95)  
  
  if(plot){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(X_comb[,1],X_comb[,2],col=X_comb$color,main="",xlim=c(xr1,xr2),ylim=c(yr1,yr2),xlab=expression(paste("Predicted"," ",NO[2]," ","(Pareto Scale)")),ylab=expression(paste("Observed"," ",NO[2]," ","(Pareto Scale)")), cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3)
    X_comb=X_comb[-a,]
    draw.circle(x = 0, y = 0, radius = R_95, border = "red", lwd=1, lty=1)
    if(Plot_Point){
      points(x = Mtx_P[a,1],y=Mtx_P[a,2],pch="*",col="blue",cex=3,lwd=2)    
    }
    if(ray){
      for(i in 1:length(Ang_T)){
        abline(a = 0, b = tan(Ang_T[i]), col="red", lty=3, lwd=0.7)  
      }
    }
    if(cone){
      abline(a = 0, b = tan(Ang_L), col="blue", lwd=3, lty=2)
      abline(a = 0, b = tan(Ang_U), col="blue", lwd=3, lty=2)
    }
  }  
  
  return(list("Ang_ecdf"=Ang_ecdf,"Mtx_P"=Mtx_P,"Mtx_top"=Mtx_top,"Total_mass"=Total_mass,"Ang_L"=Ang_L,"Ang_U"=Ang_U,"Coverage"=Coverage))
}




