#' Polar joint prediction region
#' 
#' @description
#' `jointRegion()` takes the normalized angular measure obtained from the completely positive decomposition `CPfactor()`
#' Thresholding of \eqn{\|(\widehat{X}_{p+1},X_{p+1})\|_2} at a high quantile, the function draw the bounds at the 0.025
#' and 0.975 empirical quantiles of the univariate distribution of angles provided by the normalized angular measure.
#' It returns a scatterplot of \eqn{\hat{X}_{p+1},X_{p+1}} over a high threshold with 
#' the 95\% joint region in polar coordinates. 
#' 
#' 
#' @param Xf Numeric vector; a vector of predictand
#' @param Xhat Numeric vector; a vector of the best linear prediction
#' @param Angular Numeric vector; a vector of angular components obtained from a CP-factorization of the prediction TPDM 
#' @param Pmass Numeric vector; a vector of point masses calculated from a CP-factorization of the prediction TPDM
#' @param Quan Numeric; a high quantile for radial components of \eqn{(\hat{X}_{p+1},X_{p+1})}
#' @param plot Logical; return a scatterplot with the joint region (default)
#' @param axisLimit Numeric; the max range of axes
#' @param dataPoint Numeric; a location for a specific large data point
#'
#' @return a list of quantities;
#'  accumulate angular point measure, a data frame of \eqn{(\hat{X}_{p+1},X_{p+1})},
#'  threshold excceedances of the data frame, total mass of the sum of angular measures,
#'   0.025 qunatile for w, 0.975 quantile for w, and coverage rate
#' @export
#'
jointRegion <- function(Xhat, Xf, Angular, Pmass, Quan=0.95, Plot=TRUE, axisLimit=40, dataPoint){
  
  ##  Angular measure from CP-factorization
  Ang_c=cbind(Angular,Pmass)
  sl=sort.list(Angular)
  Ang_sort=Ang_c[sl,]
  Cum_s=cumsum(Ang_sort[,2])
  Ang_H=cbind(Ang_sort,Cum_s)
  
  ##  Total mass
  Tmass=sum(Pmass)
  
  ##  Find lower and upper quantiles of angular
  Lquan=(1-Quan)/2
  Uquan=1-Lquan
  wQuan=approx(y=Ang_H[,1], x=Ang_H[,3], xout=c(Tmass*Lquan,Tmass*Uquan), method = "linear", rule=2)$y
  
  ##  Transforming (Xhat, Xf) to polar coordinates
  XY <- cbind(Xhat,Xf)
  Rad <- sqrt(apply(XY^2, 1, sum))
  W <- XY/Rad
  Keep <- Rad > quantile(Rad, Quan)
  W_top <- W[Keep,]
  XY_top <- XY[Keep,]
  
  interiorRegion=W_top[W_top[,1]>=wQuan[1] & W_top[,1]<=wQuan[2], 1]
  coverage=length(interiorRegion)/dim(W_top)[1]
  
  ##  Plot of (Xhat, Xf) with a joint polar region
  X_comb=data.frame(XY,Rad)
  X_comb$color=as.character(cut(X_comb[,3],breaks = c(0,as.numeric(quantile(X_comb[,3],Quan)),Inf),labels = c("grey","black"),right = FALSE))
  R_thres=quantile(Rad,probs = Quan)  
  
  if(Plot){
    dev.new()
    #pdf("/home/leej40/Documents/extlinear/simJointRegion.pdf",6,6)
    par(mar=c(5.1,5.1,2,2))
    plot(X_comb[,1],X_comb[,2],col=X_comb$color,main="",xaxs='i',yaxs='i',
         xlim=c(0,axisLimit),ylim=c(0,axisLimit),
         xlab=expression(paste("Predicted"," ",NO[2]," ","(Pareto Scale)")),
         ylab=expression(paste("Observed"," ",NO[2]," ","(Pareto Scale)")),
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)
    X_comb=X_comb[-dataPoint,]
    points(x = XY[dataPoint,1],y=XY[dataPoint,2],pch="*",col="blue",cex=4,lwd=2)
    draw.circle(x = 0, y = 0, radius = R_thres, border = "red", lwd=1, lty=1)
    abline(a = 0, b = sqrt(1-wQuan[1]^2)/wQuan[1], col="blue", lwd=2, lty=2)
    abline(a = 0, b = sqrt(1-wQuan[2]^2)/wQuan[2], col="blue", lwd=2, lty=2)
    #dev.off()
    }
  
  return(list("Ang_H"=Ang_H,"XY"=XY,"XY_top"=XY_top,"Tmass"=Tmass,"wQuan"=wQuan,"coverage"=coverage))
}




