#Simulation scheme

simulate <- function(n=500,time_base=c(0,0.5,1,1.5,2,2.5),
                     X_dist=c('binom'),
                     alpha0=c(log(2)),
                     alpha1=0,
                     Y_dist=c('poi','poi','bin','bin','normal','normal'),
                     beta0 = matrix(c(0.7,0.7,
                                      3,3,
                                      -0.5,-0.5,
                                      0.5,0.5,
                                      20,20,
                                      5,5),ncol=2,byrow = T),
                     beta1 = matrix(c(1,0,
                                      0,-1,
                                      1,0,
                                      0,-1,
                                      0,0,
                                      0,0.25*20),ncol=2,byrow = T),
                     phi=matrix(c(1,1,
                                  1,1,
                                  1,1,
                                  1,1,
                                  25,25,
                                  25,25),ncol=2,byrow=T),
                     gamma=matrix(rep(0.3,12),ncol=2)
){
  
  require(mvtnorm)
  
  num_class = length(alpha0) + 1
  # Generate latent class
  if (X_dist == 'binom'){
    Xcov = rbinom(n,size=1,p=0.5)
  }
  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (i in 2:num_class){
    p[,i] = exp(alpha0[i-1]+alpha1[i-1]*Xcov)
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)
  
  vecz <- matrix(0,nrow=n,ncol=num_class)
  for (i in 1:n){
    vecz[i,] = t(rmultinom(1,1,p[i,]))
  }
  z <- apply(vecz,1,which.max)
  
  # simulate each cluster (time and id)
  
  # the first subject
  maxsize <- length(time_base)
  clustersize <- maxsize  # can use sample(3:maxsize,1) to
  # simulate subjects with different observations
  time <- time_base[1:clustersize]
  id <- rep(1,clustersize)
  latent <- rep(z[1],clustersize)
  baselinecov <- rep(Xcov[1],clustersize)
  num_obs <- 1:clustersize
  
  # the second to nth subject
  for (i in 2:n){
    clustersize <- maxsize # sample(3:maxsize,1)
    id <- c(id,rep(i,clustersize))
    time <- c(time,time_base[1:clustersize])
    latent <- c(latent,rep(z[i],clustersize))
    baselinecov <- c(baselinecov,rep(Xcov[i],clustersize))
    num_obs <- c(num_obs,1:clustersize)
  }
  
  num_feature = length(Y_dist)
  y <- matrix(0,nrow=length(id),ncol=num_feature)
  
  for (j in 1:num_feature){
    
    # define link functions
    
    if (Y_dist[j] == 'normal'){
      mulink <- function(x) x
    }else if (Y_dist[j] == 'poi'){
      mulink <- function(x) exp(x)
      varlink <- function(mu) mu
    }else if (Y_dist[j] == 'bin'){
      mulink <- function(x) exp(x)/(1+exp(x))
    }
    
    # simulate longitudinal markers
    
    for (i in 1:n){
      mu <- mulink(t(beta0[j,z[i]] + beta1[j,z[i]]%*%t(time[id==i])))
      if(Y_dist[j] == 'normal'){
        
        # obtain the correlation matrix v with AR1 structure
        v <- phi[j,z[i]]*diag(1,length(mu)) %*% (gamma[j,z[i]]^as.matrix(dist(num_obs[id==i]))) %*% diag(1,length(mu))
        
        # simulate by rmvnorm
        y[id==i,j] = rmvnorm(1,mean=mu,sigma=v)
      }else if(Y_dist[j] == 'poi'){
        
        # obtain the correlation matrix v with AR1 structure
        v <- phi[j,z[i]]*sqrt(diag(varlink(as.vector(mu)))) %*% gamma[j,z[i]]^as.matrix(dist(num_obs[id==i])) %*% t(sqrt(diag(as.vector(varlink(mu)))))
        
        # simulate by poisimu function defined in file binsimu.r
        y[id==i,j] = poisimu2(mu,v)
      }else if(Y_dist[j] == 'bin'){
        
        # simulate by binsimu function defined in file binsimu.r
        y[id==i,j] = binsimu(mu,gamma[j,z[i]])
      }
      #cat(i,'\n')
    }
  }
  
  output <- data.frame(id=id,time=time,num_obs = num_obs,latent=latent,baselinecov=baselinecov,y=y)
  return(output)
}

