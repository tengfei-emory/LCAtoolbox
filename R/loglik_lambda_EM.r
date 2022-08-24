#score <- loglik(alpha,zeta,x,delta,t,num_class,recur_res)

loglik_EM_r <- function(alpha,zeta,x,delta,t,tevent,num_class,lambda_res,tau){
  
  require(Matrix)
  
  n <- length(t)
  d <- lambda_res$d
  d1_zeta <- lambda_res$d1_zeta
  d2_z <- lambda_res$d2_z
  
  l <- 0
  l_a <- 0
  l_z <- 0
  I_a <- 0
  I_z <- 0
  S <- 0
  
  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (c in 2:num_class){
    p[,c] = exp(cbind(1,x) %*% alpha[,c-1])
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)
  
  for (i in 1:length(tevent)){
    
    # find observations lies between previous event and the current event
    if (i < length(tevent)){
      idx <- which(t >= tevent[i] & t < tevent[i+1])
    }else{
      idx <- which(t >= tevent[i])
    }
    
    for (j in idx){
      
      xi <- c(1,x[j,])
      S_j <- 0
      
      for(c in 1:num_class){
        
        #initialize z_{il}
        zi <- c(rep(0,num_class-1), x[j,], rep(0,length(x[j,])*(num_class-1)))
        
        #update z by class
        if (c > 1){
          zi[c-1] = 1
          zi[(num_class+(c-1)*length(x[j,])):(num_class+c*length(x[j,])-1)] <- x[j,]
        }
        
        #####
        
        # update the estimating function / Hessian for alpha part
        qalpha <- rep(0,(num_class-1)*(ncol(x)+1))
        for (k in 2:num_class){
          qalpha[((k-1)*(ncol(x)+1)-1):((k-1)*(ncol(x)+1))] = qalpha[((k-1)*(ncol(x)+1)-1):((k-1)*(ncol(x)+1))] - p[j,k]*c(1,x[j,])
        }
        if (c > 1){
          qalpha[((c-1)*(ncol(x)+1)-1):((c-1)*(ncol(x)+1))] = qalpha[((c-1)*(ncol(x)+1)-1):((c-1)*(ncol(x)+1))] + c(1,x[j,])
        }
        
        Dalpha <- outer(p[j,-1],p[j,-1])
        diag(Dalpha) <- p[j,-1]*(1-p[j,-1])
        Dalpha <- kronecker(-Dalpha,outer(xi,xi))
        
        
        # update the estimating function/Hessian for zeta part
        qzeta0 <- zi*(delta[j] - as.numeric(exp(zi %*% zeta))*sum(d[1:i]))
        
        if (i == 1){
          sums1 = d1_zeta[1:i,]
        }else{
          sums1 = colSums(d1_zeta[1:i,])
        }
        
        qzeta = qzeta0 - as.numeric(exp(zi %*% zeta))*sums1
        qzeta = qzeta + delta[j]*d1_zeta[i,]/d[i]
        
        
        if (i == 1){
          sums2 = d2_z[,,1:i]
        }else{
          sums2 = rowSums(d2_z[,,1:i],dims=2)
        }
        
        Dzeta <-  - outer(zi,zi)*as.numeric(exp(zi %*% zeta))*sum(d[1:i]) -
          as.numeric(exp(zi %*% zeta))*outer(zi,sums1) -
          as.numeric(exp(zi %*% zeta))*outer(sums1,zi) -
          as.numeric(exp(zi %*% zeta))*sums2 -
          delta[j]*(outer(d1_zeta[i,],d1_zeta[i,]))/(d[i]^2) +
          delta[j]*d2_z[,,i]/d[i]
        
        #####
        
        #update contribution to the score function and information matrix
        
        l = l + tau[j,c]*((log(d[i]) + zi %*% zeta)*delta[j] - as.numeric(exp(zi %*% zeta))*sum(d[1:i])) + tau[j,c]*log(p[j,c])
        
        l_a = l_a + tau[j,c]*qalpha
        l_z = l_z + tau[j,c]*qzeta
        
        I_a = I_a + tau[j,c]*Dalpha
        I_z = I_z + tau[j,c]*Dzeta
        
        S_j = S_j + tau[j,c]*c(qalpha,qzeta)
      }
      
      S = S + outer(S_j,S_j)
      
    }
  }
  
  score <- c(l_a,l_z)
  
  infomat <- bdiag(I_a,I_z)
  
  output <- list()
  output$loglik = l
  output$score = score
  output$I = infomat
  output$S = S
  
  return(output)
}
