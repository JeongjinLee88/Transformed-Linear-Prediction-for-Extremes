#' Kernel density estimation in terms of angular components with boundary bias adjusted
#' 
#' @description
#' Returns a kernel density estimation to approximate an angular density function.
#' Angular components are obtained from a CP-factorization. As the angular density is bounded on \eqn{[0,1]},
#' we adjust the boundary bias using a probit transformation.
#' 
#' 
#' @param Ang Numeric vector of angular components obtained from a CP-factorization
#' @param Pmass Numeric vector of point masses from a multiple \eqn{p\times q} matrix B
#' @param bw Logical; specify a bandwidth (default:FALSE)
#' @param h Numeric; a bandwidth
#' @param plot Logical; returns a plot of the kernel density
#'
#' @return Numeric vector of kernel density estimates
#' @export
KDE_w <- function(Ang, Pmass, bw=FALSE, h=0.1, Plot=TRUE){
  
  wseq=seq(0,1,len=5000)
  
  ##  Transformed kde for h(w) at wseq
  ##  probitlink: use 'bvalue' to prevent numerical instability at close to 0 or 1.
  suppressWarnings(
  if(bw){
    kde_trans=ks::kde(x = probitlink(Ang,bvalue = .Machine$double.eps), w = Pmass, eval.points = probitlink(wseq,bvalue = .Machine$double.eps), h=h)
  }else{
    kde_trans=ks::kde(x = probitlink(Ang,bvalue = .Machine$double.eps), w = Pmass, eval.points = probitlink(wseq,bvalue = .Machine$double.eps))  
  }
  )
  
  ##  Adjust boundary bias
  kde_trans$eval.points=probitlink(kde_trans$eval.points, bvalue = .Machine$double.eps, inverse = T)
  kde_trans$estimate=(kde_trans$estimate)*(1/dnorm(kde_trans$eval.points))
  
  ##  Plot a kernel density in terms of w_L2
  if(Plot){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(kde_trans,xlim=c(0,1),main="",xlab="w",ylab="",
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    points(Ang,Pmass)
  }
  
  kde_w=kde_trans$estimate
  return(kde_w)
}
