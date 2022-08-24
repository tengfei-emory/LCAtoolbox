pscore <- function(para,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM){
  
  require(Matrix)
  
  n <- nrow(x0)
  
  l_a <- 0
  l_z <- 0
  alpha <- para[1:((num_class-1)*(ncol(x0)+1))]
  zeta <- para[(length(alpha)+1):length(para)]
  alpha <- matrix(alpha,nrow=ncol(x0)+1,ncol=num_class-1)
  
  for (i in 1:100){
    temp <- lambda_res$d
    taures = postweight(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d)
    tau = taures$tau
    lambda_res <- lambda_EM(zeta,x0,t,tevent,num_class,tau)
    # res <- loglik_EM(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d,lambda_res$d1,lambda_res$d2,tau1)
    if(sum(abs(lambda_res$d - temp)) < tolEM) break
    # if(i == 1){
    #   l0 = res$obsloglik
    # }else if(i == 2){
    #   diffl = res$obsloglik - l0
    #   l0 = res$obsloglik
    # }else if(i == 3){
    #   a = (res$obsloglik - l0)/(diffl+1e-5)
    #   diffl = res$obsloglik - l0
    #   lA = l0 + diffl/(1-a)
    #   l0 = res$obsloglik
    # }else if(i >= 4){
    #   a = (res$obsloglik - l0)/(diffl+1e-5)
    #   diffl = res$obsloglik - l0
    #   difflA = l0 + diffl/(1-a) - lA
    #   lA = l0 + diffl/(1-a)
    #   l0 = res$obsloglik
    #   if (abs(difflA) < tolEM) break
    # }
  }
  
  d = lambda_res$d
  tau = taures$tau
  p = taures$p
  d1_zeta = lambda_res$d1
  # 
  # lambda_res = lambda_EM(zeta,x0,t,tevent,num_class,tau)
  # d = lambda_res$d
  # d1_zeta = lambda_res$d1
  # 
  # if(N == 0){
  #   taures = postweight(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d)
  #   tau = taures$tau
  #   p = taures$p
  # }else{
  #   for (j in 1:N){
  #     taures = postweight(alpha,zeta,x0,delta,t,tevent,num_class,lambda_res$d)
  #     tau = taures$tau
  #     p = taures$p
  #     lambda_res = lambda_EM(zeta,x0,t,tevent,num_class,tau)
  #     d = lambda_res$d
  #     d1_zeta = lambda_res$d1
  #   }
  # }
  
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
      
      xi <- c(1,x0[j,])
      
      for(c in 1:num_class){
        
        #initialize z_{il}
        zi <- c(rep(0,num_class-1), x0[j,], rep(0,length(x0[j,])*(num_class-1)))
        
        #update z by class
        if (c > 1){
          zi[c-1] = 1
          zi[(num_class+(c-1)*length(x0[j,])):(num_class+c*length(x0[j,])-1)] <- x0[j,]
        }
        
        #####
        
        # update the estimating function / Hessian for alpha part
        qalpha <- rep(0,(num_class-1)*(ncol(x0)+1))
        for (k in 2:num_class){
          qalpha[((k-1)*(ncol(x0)+1)-ncol(x0)):((k-1)*(ncol(x0)+1))] = qalpha[((k-1)*(ncol(x0)+1)-ncol(x0)):((k-1)*(ncol(x0)+1))] - p[j,k]*c(1,x0[j,])
        }
        if (c > 1){
          qalpha[((c-1)*(ncol(x0)+1)-ncol(x0)):((c-1)*(ncol(x0)+1))] = qalpha[((c-1)*(ncol(x0)+1)-ncol(x0)):((c-1)*(ncol(x0)+1))] + c(1,x0[j,])
        }
        
        # update the estimating function/Hessian for zeta part
      
        qzeta = delta[j]*(zi + d1_zeta[i,]/d[i])
        
    
        #####
        
        #update contribution to the score function and information matrix
        
        l_a = l_a + tau[j,c]*qalpha
        l_z = l_z + tau[j,c]*qzeta
        
      }
    }
  }
  
  score = c(l_a,l_z)
  return(score)
}
