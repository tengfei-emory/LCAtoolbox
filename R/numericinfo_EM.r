numericinfo_EM <- function(l,alpha,zeta,x,eps,delta,t,tevent,num_class,tau){
  I <- matrix(0,nrow=length(alpha)+length(zeta),ncol=length(alpha)+length(zeta))

  for(i in 1:nrow(I)){
    for(j in i:ncol(I)){
      par <- c(as.vector(alpha),zeta)
      
      par[i] = par[i] + eps
      par[j] = par[j] + eps
      a <- matrix(par[1:length(alpha)],nrow=nrow(alpha),ncol=ncol(alpha))
      z <- par[-(1:length(alpha))]
      
      #tau1[tau1 < 1e-8] = 1e-8
      lambda_res <- lambda_EM(z,x,t,tevent,num_class,tau)
      tau1 <- postweight(a,z,x,delta,t,tevent,num_class,lambda_res$d)
      lambda_res <- lambda_EM(z,x,t,tevent,num_class,tau1$tau)
      res <- loglik_EM(a,z,x,delta,t,tevent,num_class,lambda_res$d,lambda_res$d1,lambda_res$d2,tau1$tau)
      l1 <- res$loglik
      
      par <- c(as.vector(alpha),zeta)
      
      par[i] = par[i] + eps
      par[j] = par[j] - eps
      a <- matrix(par[1:length(alpha)],nrow=nrow(alpha),ncol=ncol(alpha))
      z <- par[-(1:length(alpha))]
      lambda_res <- lambda_EM(z,x,t,tevent,num_class,tau)
      tau1 <- postweight(a,z,x,delta,t,tevent,num_class,lambda_res$d)
      lambda_res <- lambda_EM(z,x,t,tevent,num_class,tau1$tau)
      res <- loglik_EM(a,z,x,delta,t,tevent,num_class,lambda_res$d,lambda_res$d1,lambda_res$d2,tau1$tau)
      l2 <- res$loglik
      
      par <- c(as.vector(alpha),zeta)
      
      par[i] = par[i] - eps
      par[j] = par[j] + eps
      a <- matrix(par[1:length(alpha)],nrow=nrow(alpha),ncol=ncol(alpha))
      z <- par[-(1:length(alpha))]
      lambda_res <- lambda_EM(z,x,t,tevent,num_class,tau)
      tau1 <- postweight(a,z,x,delta,t,tevent,num_class,lambda_res$d)
      lambda_res <- lambda_EM(z,x,t,tevent,num_class,tau1$tau)
      res <- loglik_EM(a,z,x,delta,t,tevent,num_class,lambda_res$d,lambda_res$d1,lambda_res$d2,tau1$tau)
      l3 <- res$loglik
      
      I[i,j] = -(l1 - l2 - l3 + l)/(eps^2)
    }
  }
  
  dI = diag(I)
  
  I = I + t(I) - diag(dI)
  
  return(I)
  
}