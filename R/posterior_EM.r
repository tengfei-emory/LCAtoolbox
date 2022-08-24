postweight_r <- function(alpha,zeta,lambda_res,x,tevent,t,delta,num_class){
  
  n <- length(delta)
  
  survlik <- matrix(1,nrow=n,ncol=num_class)
  tau <- matrix(0.5,nrow=n,ncol=num_class)
  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (c in 2:num_class){
    p[,c] = exp(cbind(1,x) %*% alpha[,c-1])
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)
  d = lambda_res$d
  
  for (i in 1:length(tevent)){
    # find observations lies between previous event and the current event
    if (i < length(tevent)){
      idx <- which(t >= tevent[i] & t < tevent[i+1])
    }else{
      idx <- which(t >= tevent[i])
    }
    for (j in idx){
      
      xi <- c(1,x[j,])
      
      for(c in 1:num_class){
        
        #initialize z_{il}
        zi <- c(rep(0,num_class-1), x[j,], rep(0,length(x[j,])*(num_class-1)))
        
        #update z by class
        if (c > 1){
          zi[c-1] = 1
          zi[(num_class+(c-1)*length(x[j,])):(num_class+c*length(x[j,])-1)] <- x[j,]
        }
        
        #update posterior probability for the jth observation
        
        survlik[j,c] <- as.numeric(exp(zi %*% zeta))^(delta[j])*exp(
          -as.numeric(exp(zi %*% zeta))*sum(d[1:i]))

      }
      
      pew <- p*survlik
      #pew[pew<1e-8] = 1e-8
      pewsum = rowSums(pew)
      tau = apply(pew,2,function(x) x/pewsum)
      #tau[tau<1e-8] = 1e-8
    }
    
  }
  
  output = list()
  output$tau = tau
  output$p = p
  return(output)
  
}