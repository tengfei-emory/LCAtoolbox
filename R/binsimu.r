binsimu <- function(mu,gamma){
  
  # simulate correlated binary variables using method by Qaqish (2003)
  
  len <- length(mu)
  y <- rep(0,len)
  y[1] <- rbinom(1,1,mu[1])
  for (i in 2:len){
    lambda <- mu[i] + gamma*(y[i-1] - mu[i-1])*sqrt((mu[i]*(1-mu[i]))/(mu[i-1]*(1-mu[i-1])))
    y[i] <- rbinom(1,1,lambda)
  }
  y
}

poisimu <- function(mu,v){
  
  # simulate correlated poisson variables using the first method by Dalthorp and Madsen (2007)
  
  muz <- log(mu^2/sqrt(diag(v)-mu+mu^2))
  covz <- matrix(0,length(muz),length(muz))
  diag(covz) <- log((diag(v)-mu)/(mu^2)+1)
  
  for (i in 1:nrow(covz)){
    for (j in 1:ncol(covz)){
      if(i != j){
        covz[i,j] = log(v[i,j]/(mu[i]*mu[j])+1)
      }
    }
  }
  
  z <- rmvnorm(n=1,mean=muz,sigma=covz)
  z <- exp(z)
  y <- rep(0,length(mu))
  for (i in 1:length(mu)){
    y[i] <- rpois(1,z[i])
  }
  y
}

poisimu2 <- function(mu,v){
  l = 0
  alpha = v
  for (i in 1:length(mu)){
    for (j in 1:length(mu)){
      alpha[i,j] = alpha[i,j]#*sqrt(mu[i]*mu[j])
    }
  }
  vp <- NULL
  G = matrix(0,ncol=ncol(v),nrow=nrow(v))
  G[alpha > 0] = 1
  
  while (max(alpha) > 0){
    l = l + 1
    #cat(l,'\n')
    
    inds = which(alpha == min(alpha[alpha>0]), arr.ind=TRUE)[1,]
    beta = min(alpha[alpha>0])
    vp <- c(vp,beta)
    
    S = which(G[,inds[1]] == 1 & G[,inds[2]] == 1)
    for (i in S){
      for (j in S){
        alpha[i,j] = alpha[i,j] - beta
      }
    }
    
    alpha[alpha<1e-8] = 0
    
    if(l == 1){
      tvec = matrix(0,nrow=nrow(G),ncol=1)
      tvec[S] = 1
      Tmat = tvec
    }else{
      tvec = matrix(0,nrow=nrow(G),ncol=1)
      tvec[S] = 1
      Tmat = cbind(Tmat,tvec)
    }
    
    G[alpha == 0] = 0
  }
  
  # idx <- which(vp < 0.25)
  
  # Gmat <- matrix(0,ncol=ncol(Tmat),nrow=ncol(Tmat))
  # Gmat[1:ncol(Tmat),] = diag(ncol(Tmat))
  # Gmat = rbind(Gmat,-diag(ncol(Tmat)))
  # Hvec <- rep(0,ncol(Tmat))
  # Hvec <- c(Hvec,-vp-1)
  # for(id in idx){
  #   Gvec <- rep(0,ncol(Tmat))
  #   Gvec[id] = 1 
  #   Gmat = rbind(Gmat,Gvec)
  #   Gmat = rbind(Gmat,-Gvec)
  #   Hvec = c(Hvec,0.5 + sqrt(0.25-vp[id]))
  #   Hvec = c(Hvec,-0.5 + sqrt(0.25-vp[id]))
  # }
  
  #mup <- ldei(Tmat,mu,Gmat,Hvec)$X
  #mup <- Solve(Tmat,mu)
  #vp <- vp[mup > 1e-8]
  #Tmat <- Tmat[,mup > 1e-8]
  #mup <- mup[mup > 1e-8]
  
  X <- rep(NA,length(vp))
  # for(k in 1:length(X)){
  #   if (vp[k] > mup[k]){
  #     X[k] <- rnbinom(1,size=mup[k]^2/(vp[k]-mup[k]),prob=mup[k]/vp[k])
  #   }else if(vp[k] < mup[k]){
  #     X[k] <- rbinom(1,size=1,prob=sqrt(mup[k]-vp[k]))+rpois(1,mup[k]-sqrt(mup[k]-vp[k]))
  #   }
  # }
  for (k in 1:length(X)){
    X[k] = rpois(1,vp[k])
  }
  
  Y = Tmat%*%X
  Y
}