obsloglik_totN <- function(para,x0,delta,t,tevent,num_class,tau,lambda_res,N){
  
  require(Matrix)
  
  n <- length(t)
  
  l_ind <- rep(0,n)
  l <- 0
  
  alpha <- para[1:((num_class-1)*(ncol(x0)+1))]
  zeta <- para[(length(alpha)+1):length(para)]
  alpha <- matrix(alpha,nrow=ncol(x0)+1,ncol=num_class-1)
  
  if(N == 0){
    taures = postweight(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d)
    tau = taures$tau
    p = taures$p
    d = lambda_res$d
  }else{
    for (j in 1:N){
      taures = postweight(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d)
      tau = taures$tau
      p = taures$p
      lambda_res = lambda_EM(zeta,x0,t,tevent,num_class,tau)
      d = lambda_res$d
    }
  }
  
  for (i in 1:length(tevent)){
    
    # find observations lies between previous event and the current event
    
    if (i == 1){
      idx <- which(t < tevent[i+1])
    }else if (i < length(tevent)){
      idx <- which(t >= tevent[i] & t < tevent[i+1])
    }else{
      idx <- which(t >= tevent[i])
    }
    
    for (j in idx){
      
      for(c in 1:num_class){
        
        #initialize z_{il}
        zi <- c(rep(0,num_class-1), x0[j,], rep(0,length(x0[j,])*(num_class-1)))
        
        #update z by class
        if (c > 1){
          zi[c-1] = 1
          zi[(num_class+(c-1)*length(x0[j,])):(num_class+c*length(x0[j,])-1)] <- x0[j,]
        }
        
        #update contribution to the score function and information matrix
        
        ds = sum(d[1:i])
        
        if(t[j] < tevent[i]){
          l_ind[j] = l_ind[j] + p[j,c]
        }else{
          l_ind[j] = l_ind[j] + p[j,c]*((d[i]*exp(zi%*%zeta))^delta[j])*exp(-ds*exp(zi%*%zeta))
        }
      }
    }
  }
  
  l = sum(log(l_ind))
  
  return(l)
}
