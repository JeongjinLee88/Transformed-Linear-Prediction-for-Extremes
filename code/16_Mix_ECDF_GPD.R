Mix_ECDF_GPD <- function(x)
{
  z <- numeric(length(x))
  below <- x <= thold
  
  shift=0.9352074
  
  ## below threshold
  zBelow <- z[below]
  u <- approxExtrap(x = CDF_c[,1],y = CDF_c[,2],xout = x[below])$y
  u_b <- u > 0 & u < 1 & !is.na(u)
  zBelow[u_b] <- 1/sqrt(1-u[u_b])-shift
  zBelow[!u_b] <- NA
  z[below]<-zBelow
    
  ##  above threshold
  zAbove <- z[!below]
  u <- 1 - rate*(1 - pgpd(x[!below], thold, scale, shape))
  u_u <- u < 1
  zAbove[u_u] <- qgpd(u[u_u], loc = 1, scale = 0.5, shape = 0.5)-shift
  #zAbove[u_u] <- ((1-u[u_u])^(-0.5)-1)+1-shift #same as above
  zAbove[!u_u] <- Inf
  z[!below] <- zAbove

  return(z)
}
