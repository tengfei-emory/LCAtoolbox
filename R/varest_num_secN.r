VarEst_num_secN <- function(para,x0,delta,t,tevent,num_class,tau,lambda_res,N){
  
  n <- length(t)
  
  h = 5*n^(-0.5)  #### tuning parameter
  
  p = length(para)
  
  E = h*diag(rep(1,p))  ### p by p
  
  Me = matrix(0,p,p)
  
  p0 = obsloglik_totN(para,x0,delta,t,tevent,num_class,tau,lambda_res,N)
  pi = rep(0,p)
  for (i in 1:p){
    paraEi = para + E[,i]
    pi[i] = obsloglik_totN(paraEi,x0,delta,t,tevent,num_class,tau,lambda_res,N)
  }
  
  for (i in 1:p){
    for (j in i:p){
      paraEij = para + E[,i] + E[,j]
      
      # p0 = obsloglik_tot(para,x0,delta,t,tevent,num_class,tau,tolEM)
      pi1 = pi[i]
      pj1 = pi[j]
      pij2 = obsloglik_totN(paraEij,x0,delta,t,tevent,num_class,tau,lambda_res,N)
      
      Me[i,j] = (p0-pi1-pj1+pij2)/(h^2) ####  score function, n by 1
      Me[j,i] = Me[i,j]
    }
  }
  
  
  #I = Me%*%t(Me)
  
  varmat = sqrt(diag(solve(-Me)))
  
  varmat
}