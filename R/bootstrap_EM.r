bootstrap_EM <- function(B=10,alpha0,zeta0,baseline,num_class,tolEM,tolNewton,maxiter,initial){
  
  b.alpha <- matrix(0,nrow=B,ncol=length(alpha0))
  b.zeta <- matrix(0,nrow=B,ncol=length(zeta0))
  n <- nrow(baseline)
  
  for(b in 1:B){
    r.id <- sample(1:n,n,replace=T)
    b.dat <- baseline[r.id,]
    
    # extract baseline data
    b.baseline <- b.dat[order(b.dat$tildet,b.dat$id,b.dat$time),]
    
    x <- as.matrix(b.baseline$baselinecov)
    delta <- b.baseline$delta
    t <- b.baseline$tildet
    tevent <- b.baseline$tildet[baseline$delta]
    
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    if (initial == 'random'){
      prop <- 25*num_class/n
      u <- runif(n,0,1)
      for (i in 1:n){
        if (u[i] < prop){
          tau0[i,ceiling(num_class*runif(1))] = 1
        }
      }
      tau0[tau0 < 1e-8] = 1e-8
    }else if(initial == 'true'){
      for (i in 1:n){
        if(baseline$latent[i] == 1){
          tau0[i,] = c(1,0)
        }else{
          tau0[i,] = c(0,1)
        }
      }
    }
    
    diff_EM = 10
    
    alpha = alpha0
    zeta = zeta0
    
    for (count in 1:maxiter){
      diff = 10
      c = 1
      while (diff > tolNewton & c < 10){
        c = c+1
        lambda_res <- Lambda_EM(alpha,zeta,x,t,tevent,num_class,tau0)
        res <- loglik_EM(alpha,zeta,x,delta,t,tevent,num_class,lambda_res,tau0)
        step <- solve(res$I+1e-4)%*%res$score
        if (sum(is.na(step)) > 0) break
        alpha = alpha - step[1:length(alpha)]
        zeta = zeta - step[(length(alpha)+1):(length(alpha)+length(zeta))]
        diff = max(abs(step))
      }      

      lambda_res <- Lambda_EM(alpha,zeta,x,t,tevent,num_class,tau0)
      postcalc <- postweight(alpha,zeta,lambda_res,x,tevent,t,delta,num_class)
      tau1 <- postcalc$tau
      tau1[tau1 < 1e-8] = 1e-8

      p <- postcalc$p
      #diff_EM = sum(abs((tau1 - tau0)^2))
      diff_EM = max(abs(c(alpha0-alpha,zeta0-zeta)))
      tau0 <- tau1
      
      if(diff_EM < tolEM) break
    }
    
    b.alpha[b,] = alpha
    b.zeta[b,] = zeta
  }
  
  ASE <- c(apply(b.alpha,2,function(x) sqrt(var(x))),apply(b.zeta,2,function(x) sqrt(var(x))))
    
  return(ASE)
}