#' simulation.lcphR: simulate a dataset for latent class survival analysis
#' @description Simulate a dataset according to the scenario (I) in the manuscript.
#' @param n sample size of the simulated dataset.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Fei, T., Hanfelt, J., & Peng, L. (2022). Latent Class Analysis with Semi-parametric Proportional Hazards Submodel for Time-to-event Data. arXiv preprint arXiv:2202.00775.
#' @examples 
#' 
#' dat <- simulation.lcphR(n=1000)
#' 
#' @export

simulation.lcphR <- function(n=1000){
  
  time_base=time_base=c(0,0.5,1,1.5,2,2.5)
  X_dist=c('binom')
  alpha0=c(log(2))
  alpha1=c(0)
  alpha2=c(0)
  zeta0=c(0,2)
  zeta1=c(-2,0)
  zeta2=c(0,2)
  lambda=0.1
  lunif=5
  exprate=0.1
  
  #require(SimCorrMix)
  require(mvtnorm)
  #require(PoisNor)
  #require(GenOrd)
  
  num_class = length(alpha0) + 1
  # Generate latent class
  if (X_dist == 'binom'){
    Xcov = rbinom(n,size=1,p=0.5)
  }
  
  Xcov2 = runif(n,0,1)
  
  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (i in 2:num_class){
    p[,i] = exp(alpha0[i-1]+alpha1[i-1]*Xcov+alpha2[i-1]*Xcov2)
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)
  
  vecz <- matrix(0,nrow=n,ncol=num_class)
  for (i in 1:n){
    vecz[i,] = t(rmultinom(1,1,p[i,]))
  }
  z <- apply(vecz,1,which.max)
  
  # simulate observed survival time
  tildet <- rep(NA,n)
  delta <- rep(NA,n)
  for(i in 1:n){
    U <- runif(1,min=0,max=1)
    #t <- -2 + 2*(1-U)^(-1/(0.5*exp(zeta0[z[i]]+zeta1[z[i]]*Xcov[i])))
    t <- log( 1-(log(1-U))/(lambda* exp(zeta0[z[i]]+zeta1[z[i]]*Xcov[i] + zeta2[z[i]]*Xcov2[i] ) ) )
    c <- min(rexp(1,rate=exprate),runif(1,min=lunif,max=lunif+1))
    tildet[i] <- min(t,c)
    delta[i] <- I(t <= c)
  }
  
  output <- data.frame(id=1:n,time=tildet,latent=z,Xcov1=Xcov,Xcov2=Xcov2,tildet=tildet,delta=delta)
  return(output)
}