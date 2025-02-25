#' Coverage rates for prediction interval
#' 
#' @description
#' Returns a coverage rate for the 95\% prediction interval, for example.
#' The approximate conditional density given a large predicted value is obtained from 
#' the function `conDensity()`.
#' 
#' @param XY \eqn{n\times 2} numeric matrix of \eqn{\hat{x}_{p+1}} and \eqn{x_{p+1}}
#' @param kde_est Numeric vector of kernel density estimates
#' @param Plot Logical; returns a plot of \eqn{\hat{x}_{p+1}} and \eqn{x_{p+1}} with the 95\% prediction interval
#'
#' @return a list of numeric value of coverage rate and a matrix of prediction intervals
#' @export
coverageRate <- function(XY, kde_est, Quan, Plot=TRUE){
  
  XY_top=XY

  ##  Assess the coverage rate
  counting=rep(0,dim(XY_top)[1])
  PI=matrix(0,nrow = dim(XY_top)[1],ncol = 2,byrow = T)
  for(i in 1:dim(XY_top)[1]){
    out=condDensity(xhatPoint=XY_top[i,1],xp1Point=XY_top[i,2],
                    kde_h = kde_est, Quan = Quan,
                    Plot_v = FALSE, Plot_h = FALSE)
    counting[i]=out$count
    PI[i,]=out$condQuan
    print(i)
  }
  CoverageRate=sum(counting)/dim(XY_top)[1]
  
  ##  Plot of (Xhat,X_p+1)
  ABC=cbind(XY_top[,1],PI)
  sl <- sort.list(ABC[,1])
  ABCSort=ABC[sl,]
  if(Plot){
    dev.new()
    #pdf("/home/leej40/Documents/extlinear/no2ParetoScale.pdf",6,6)
    #pdf("/home/leej40/Documents/extlinear/simCondIntvls.pdf",6,6)
    par(mar=c(5.1,5.1,2,2))
    plot(XY_top[,1], XY_top[,2], xlim=c(0,70), ylim=c(0,70),
         main="", xlab="", ylab=expression(X[p+1]),
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3)
    title(xlab=expression(hat(X)[p+1]),line=3.5, cex.lab=1.5)
    lines(ABCSort[,1],ABCSort[,2], lty = 2, cex=2, lwd=2, col="blue")
    lines(ABCSort[,1],ABCSort[,3], lty = 2, cex=2, lwd=2, col="blue")
    #points(x = XY_top[,1],y=PI[,1],pch="-", col="blue",cex=2,lwd=2)
    #points(x = XY_top[,1],y=PI[,2],pch="-", col="blue",cex=2,lwd=2)
    #dev.off()
  }
  
  return(list("CoverageRate"=CoverageRate,"PI"=PI))
}
