Ang_kde <- function(Ang_T, Ang_mass, bw=T, h=0.1, plot_theta=T, plot_w=T, y_lim){
  
  ##  Transform angles to angular components
  w_L2=cos(Ang_T)  
  
  ##  Kernel density with a bandwidth 
  if(bw){
    kde=ks::kde(x = probitlink(Ang_T/(pi/2)), w = Ang_mass, h=h)
  }
  else{
    kde=ks::kde(x = probitlink(Ang_T/(pi/2)), w = Ang_mass)      
  }
  
  kde_trans=kde
  kde_trans$eval.points=probitlink(kde$eval.points,inverse = T)*(pi/2)
  xx=probitlink(kde_trans$eval.points/(pi/2))
  
  ##  Adjust boundary bias
  kde_trans$estimate=(kde_trans$estimate)*(1/dnorm(xx))/(pi/2)
  
  ##  Plot a kernel density in terms of angles
  if(plot_theta){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    #plot(kde_trans,xlim=c(0,pi/2),ylim=c(0,max(kde_trans$estimate)*1.5),main="The weighted kernel angular density",xlab=expression(theta), cex.main=1.3, cex.lab=1.3, cex.axis=1.3)
    plot(kde_trans,xlim=c(0,pi/2),ylim=c(0,1.5),main="",xlab=expression(theta),ylab="", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    points(Ang_T,Ang_mass)
  }
  
  ##  Evaluate h(w) in terms of angular components
  theta=seq(0.0001, pi/2-0.0001, length.out=10000)
  
  w1 <- cos(theta)
  if(bw){
    kde_w=ks::kde(x = probitlink(w_L2), w = Ang_mass, eval.points = probitlink(w1), h=h)
  }
  else{
    kde_w=ks::kde(x = probitlink(w_L2), w = Ang_mass, eval.points = probitlink(w1))  
  }
  
  kde_trans_w=kde_w
  kde_trans_w$eval.points=probitlink(kde_w$eval.points,inverse = T)
  xx_w=probitlink(kde_trans_w$eval.points)
  
  ##  Adjust boundary bias
  kde_trans_w$estimate=(kde_trans_w$estimate)*(1/dnorm(xx_w))
  #kde_trans_w$estimate=(kde_trans_w$estimate)*(1/dnorm(xx_w))*(kde_trans_w$eval.points*2)
  
  ##  Plot a kernel density in terms of w_L2
  if(plot_w){
    dev.new()
    par(mar=c(5.1,5.1,2,2))
    plot(kde_trans_w,xlim=c(0,1),ylim=c(0,y_lim),main="",xlab="W",ylab="", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    points(w_L2,Ang_mass)
  }
  
  return(list("kde_trans"=kde_trans$estimate,"kde_trans_w"=kde_trans_w$estimate))
}
