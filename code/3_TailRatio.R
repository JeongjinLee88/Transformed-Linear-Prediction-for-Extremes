TailRatioPlot <- function(n, Ncol, A, plot_x=T, plot_preimage=T){
  
  ##  Simulate a bivariate Pareto dist
  U <- runif(n*Ncol)
  Z <- matrix(1/sqrt(1-U),nrow=n,ncol=Ncol)
  
  ##  A TPDM (known)
  Azero=pmax(A,0)
  TPDM=Azero%*%t(Azero)
  
  ##  Calculate a vector 'X' and its preimage 'Y'
  X <- Amul(A, t(Z))
  Y <- A %*% InvT(x = t(Z))
  if(plot_x){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(t(X), xlab = expression(X[1]), ylab = expression(X[2])
         , cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  if(plot_preimage){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(t(Y), xlab = expression(Y[1]), ylab = expression(Y[2])
         , cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  return(list("X"=t(X),"Y"=t(Y),"TPDM"=TPDM))
}




