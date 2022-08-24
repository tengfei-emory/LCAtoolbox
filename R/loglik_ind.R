loglik_ind <- function(para,x0,delta,t,tevent,num_class,tau,j){
  
  require(Matrix)
  
  n <- nrow(x0)
  
  l <- 0
  
  alpha <- para[1:((num_class-1)*(ncol(x0)+1))]
  zeta <- para[(length(alpha)+1):length(para)]
  alpha <- matrix(alpha,nrow=ncol(x0)+1,ncol=num_class-1)
  lambda_res = lambda_EM(zeta,x0,t,tevent,num_class,tau)
  d = lambda_res$d
  taures = postweight(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d)
  tau = taures$tau
  p = taures$p
  
  # p <- matrix(0,ncol=num_class,nrow=n)
  # p[,1] = 1
  # for (c in 2:num_class){
  #   p[,c] = exp(cbind(1,x0) %*% alpha[,c-1])
  # }
  # psum = rowSums(p)
  # p = apply(p,2,function(x0) x0/psum)
  
  i = sum(delta[1:j])
  
  for(c in 1:num_class){
    
    #initialize z_{il}
    zi <- c(rep(0,num_class-1), x0[j,], rep(0,length(x0[j,])*(num_class-1)))
    
    #update z by class
    if (c > 1){
      zi[c-1] = 1
      zi[(num_class+(c-1)*length(x0[j,])):(num_class+c*length(x0[j,])-1)] <- x0[j,]
    }
    
    #update contribution to the score function and information matrix
    
    l = l + tau[j,c]*log(p[j,c])
    
    if(i > 0){
      l = l + tau[j,c]*((log(d[i]) + zi %*% zeta)*delta[j])
    }
  }
  
  return(l)
}
