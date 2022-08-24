VarEst_numN <- function(para,x0,delta,t,tevent,num_class,tau,lambda_res,N){
  
  n <- length(t)
  
  h = 5*n^(-0.5)  #### tuning parameter
  
  p = length(para)
  
  E = h*diag(rep(1,p))  ### p by p
  
  Me = matrix(0,(p),n)
  
  i = 1
  
  for (i in 1:p){
    
    paraEi1 = para + E[,i]
    paraEi2 = para - E[,i]
    
    p1 = obsloglikN(paraEi1,x0,delta,t,tevent,num_class,tau,lambda_res,N)
    
    p2 = obsloglikN(paraEi2,x0,delta,t,tevent,num_class,tau,lambda_res,N)
    
    Me[i,] = (p1-p2)/(2*h) ####  score function, n by 1
    
  }
  
  
  I = Me%*%t(Me)
  
  varmat = sqrt(diag(solve(I)))
  
  varmat
}