# d1 = 1/10000
# t = baseline$tildet
# tevent = baseline$tildet[baseline$delta==1]
# x = as.matrix(baseline$baselinecov)
# zeta = c(2,-1,2)
# alpha = matrix(c(0.6,0),nrow=2,ncol=1)
# delta = baseline$delta
# num_class = 2
# recur_res <- recursive(d1,alpha,zeta,x,delta,t,tevent,num_class)

Lambda_EM <- function(alpha,zeta,x,t,tevent,num_class,tau){
  
  n <- length(t)
  
  # initializing first-order partial d and second-order partial d
  
  d <- rep(0,length(tevent))
  
  d1_zeta <- matrix(0,nrow=length(tevent),ncol=length(zeta))
  
  d2_z <- array(0,dim=c(ncol(d1_zeta),ncol(d1_zeta),length(tevent)))
  
  dall <- rep(0,n)
  d1all <- matrix(0,nrow=n,ncol=length(zeta))
  d2all <- array(0,dim=c(ncol(d1_zeta),ncol(d1_zeta),n))
  
  for (j in 1:n){
    
    for (c in 1:num_class){
      
      #initialize z_{il}
      zi <- c(rep(0,num_class-1), x[j,], rep(0,length(x[j,])*(num_class-1)))
      
      #update z by class
      if (c > 1){
        zi[c-1] = 1
        zi[(num_class+(c-1)*length(x[j,])):(num_class+c*length(x[j,])-1)] <- x[j,]
      }
      
      #update ith jump 
      
      dall[j] <- dall[j] + tau[j,c]*as.numeric(exp(zi %*% zeta))
      
      d1all[j,] <- d1all[j,] + tau[j,c]*zi*as.numeric(exp(zi %*% zeta))
      
      d2all[,,j] <- d2all[,,j] + tau[j,c]*outer(zi,zi)*as.numeric(exp(zi %*% zeta)) 
      
    }
    
  }
  
  for (i in 1:length(tevent)){
    
    # find observations lies between previous event and the current event
    idx <- which(t >= tevent[i])
    
    if(length(idx) == 1){
      d2_z[,,i] <-  d2all[,,idx]/(dall[idx]^2) - 2*outer(d1all[idx,],d1all[idx,])/(dall[idx]^3)
      d1_zeta[i,] <- d1all[idx,]/(dall[idx]^2)
      d[i] <- 1/dall[idx]
    }else{
      d2_z[,,i] <-  rowSums(d2all[,,idx],dims=2)/(sum(dall[idx])^2) - 2*outer(colSums(d1all[idx,]),colSums(d1all[idx,]))/(sum(dall[idx])^3)
      d1_zeta[i,] <- colSums(d1all[idx,])/(sum(dall[idx])^2)
      d[i] <- 1/sum(dall[idx])
    }
    
    
  }
  
  output <- list()
  output$d = d
  output$d1_zeta = d1_zeta
  output$d2_z = d2_z
  
  return(output)
}
