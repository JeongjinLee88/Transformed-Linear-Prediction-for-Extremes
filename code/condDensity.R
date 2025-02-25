#' Approximate conditional angular density given large predicted values
#' 
#' @description
#' Return an approximate conditional angular density given a large value \eqn{\hat{X}_{p+1}=\hat{x}_{p+1}}.
#' The function takes the kernel density estimate of angular masses obtained from 
#' the completely positive decomposition of the \eqn{2\times 2} prediction TPDM estimate
#' 
#' @param xhatPoint Numeric; a large prediction
#' @param xp1Point Numeric; the corresponding \eqn{x_{p+1}} value
#' @param kde_h Numeric vector; kernel estimates based on angular components and their point masses via CP-factorization
#' @param Quan Numeric; a high quantile
#' @param xlim_upper Numeric; the upper limit of the x-axis
#' @param ylim_upper Numeric; the upper limit of the y-axis
#' @param Plot_v Logical; returns the approximate conditional density given \eqn{\hat{x}_{p+1}} (default:TRUE)
#' @param Plot_h Logical; returns the approximate conditional density given \eqn{\hat{x}_{p+1}} with the 95\% conditional interval horizontally
#'
#' @return a list of vectors; approximated conditional density
#' @export
condDensity <- function(xhatPoint, xp1Point, kde_h, Quan, xlim_upper, ylim_upper, Plot_v = TRUE, Plot_h = TRUE){
  
  ##  Calculate the approximate conditional density for a given large Xhat
  ##  Define angular components 
  wseq=seq(1,1e-10,len=length(kde_h))
  ##  Radial components and corresponding x_p+1
  rad <- xhatPoint/wseq
  #xp1 <- rad*sqrt(1-wseq^2)
  xp1 <- rad*(1-wseq)
  
  #xp1<-xhatPoint*sqrt((1-wseq^2)/(wseq^2))
  #rad<- sqrt(xhatPoint^2+xp1^2)
  
  ##  'denApproxNorm' returns an approximated angular density from CPD
  denApproxNoNorm=2*(rad^(-5))*rev(kde_h)*xp1
  denApproxNoNorm[is.nan(denApproxNoNorm)] <- 0
  
  xp1Temp <- xp1
  #xp1Temp[length(xp1Temp)] <- 2 * xp1Temp[length(xp1Temp) - 1]
  delta <- xp1Temp[2:(length(xp1Temp))] - xp1Temp[1:(length(xp1Temp) - 1)]
  piece1 <- delta * denApproxNoNorm[1:(length(denApproxNoNorm) - 1)]
  piece2 <- delta * denApproxNoNorm[2:(length(denApproxNoNorm))]
  traps <- 1/2*(piece1 + piece2)
  normalizer <- sum(traps)
  cumTraps <- c(0, cumsum(traps))/normalizer # CDF
  denApproxNorm <- denApproxNoNorm/normalizer
  
  
  ##  Sort the CDF
  Con_comb=cbind(xp1,cumTraps)
  sl=sort.list(xp1)
  Con_sort=Con_comb[sl,]
  ##  Interpolate the empirical CDF to find the 95% conditional interval
  Lquan=(1-Quan)/2;Uquan=1-Lquan
  condQuan=approx(x = Con_sort[,2], y = Con_sort[,1], xout = c(Lquan,Uquan), method = "linear",rule = 2,ties = "mean")$y
  #condQuan=approxExtrap(x = Con_sort[,2], y = Con_sort[,1], xout = c(Lquan,Uquan), method = "linear",rule = 2,ties = "mean")$y
  
  ##  Counting
  count=ifelse(xp1Point>=condQuan[1] & xp1Point<=condQuan[2],1,0)
  
  ##  Plot of the approximate conditional density with the 95% CI

  if(Plot_h){
    #dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(xp1, denApproxNorm, type = 'l', 
         ylim = range(denApproxNorm), xlim=c(0,50),
         xlab=expression(X[p+1]~"|"~hat(X)[p+1]),
         ylab = "Conditional density",
         main=expression(paste("The"," ","approximated"," ","conditional"," ","density"," ","when"," ",hat(X)[p+1]," ","is"," ","large")))
  }
  if(Plot_v){
    #dev.new()
    #pdf("/home/leej40/Documents/extlinear/simCondIntvlXHat27.pdf",6,6)
    par(mar=c(5.1,5.1,2,2))
    plot(xhatPoint+denApproxNorm, xp1, type='l', lwd=2,
         xlim=c(xhatPoint-0.01,xlim_upper),ylim=c(0,ylim_upper),
         xlab="Predicted values", ylab="Observed values",
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)
    abline(v = xhatPoint,lwd=2)
    points(x = xhatPoint, y=xp1Point, pch="*", col="blue", cex=3, lwd=2)
    ##  95% conditional interval  
    points(x = xhatPoint,y=condQuan[1],pch="-", col="red",cex=3,lwd=2)
    points(x = xhatPoint,y=condQuan[2],pch="-", col="red",cex=3,lwd=2)
    #dev.off()
  }
  
  return(list("denApproxNoNorm"=denApproxNoNorm, "denApproxNorm"=denApproxNorm, "xp1"=xp1, "cumTraps"=cumTraps, "condQuan"=condQuan, "count"=count))
}

