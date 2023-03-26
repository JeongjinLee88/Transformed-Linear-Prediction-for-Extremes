CondDens <- function(z, h_w2_CPD, plot = T){
  
  ##  The approximate conditional density when Xhat is large
  ##  z is a given large value (Xhat)
  z1 <- z[1]
  
  ##  Define angular components 
  #theta=seq(0.0000001, pi/2-0.0000001, length.out = 501)
  theta=seq(0.0001, pi/2-0.0001, length.out=10000)
  w1 <- cos(theta)
  w2 <- sin(theta)
  w <- cbind(w1, w2)
  
  ##  Gives corresponding values of angular density from CPD
  sd_CPD <- h_w2_CPD
  ##  Gives corresponding radial components
  rad_CPD <- z1/w1
  ##  Gives corresponding values of z2
  z2_CPD <- w2*rad_CPD
  
  ##  The resulting 'densVal1' returns an approximated angular density from CPD 
  unDensVal=(rad_CPD^(-5))*sd_CPD*z2_CPD
  unDensVal[is.nan(unDensVal)] <- 0
  
  z2Temp_CPD <- z2_CPD
  z2Temp_CPD[length(z2Temp_CPD)] <- 2 * z2Temp_CPD[length(z2Temp_CPD) - 1]
  deltaz3 <- z2Temp_CPD[2:(length(z2Temp_CPD))] - z2Temp_CPD[1:(length(z2Temp_CPD) - 1)] 
  piece1 <- deltaz3 * unDensVal[1:(length(unDensVal) - 1)]
  piece2 <- deltaz3 * unDensVal[2:(length(unDensVal))]
  traps <- 1/2*(piece1 + piece2)
  normalizer <- sum(traps)
  cumTraps <- c(0, cumsum(traps))/normalizer
  densVal <- unDensVal/normalizer
  
  if(plot){
    dev.new()
    plot(z2_CPD, densVal, type = 'l', xlim = c(0,2*sum(z)), ylab = "Conditional density", ylim = range(densVal), main=expression(paste("The"," ","approximated"," ","conditional"," ","density"," ","when"," ",hat(X)[p+1]," ","is"," ","large"))
         ,xlab=expression(X[p+1]~"|"~hat(X)[p+1]))
    legend("topright",legend=c("Approx Cond dens from CPD"),lty=c(3,1))
  }
  
  return(list("unDensVal"=unDensVal, "densVal"=densVal, "z2_CPD"=z2_CPD, "cumTraps"=cumTraps))
}

