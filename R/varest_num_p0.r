VarEst_num_p0 <- function(para,x0,delta,t,tevent,num_class,tau,tolEM){
  
  n <- length(t)
  
  h = 5*n^(-0.5)  #### tuning parameter
  
  p = length(para)
  
  E = h*diag(rep(1,p))  ### p by p
  
  Me = matrix(0,(p),n)
  
  i = 1
  
  p0 = obsloglik(para,x0,delta,t,tevent,num_class,tau,tolEM)
  
  for (i in 1:p){
    
    paraEi1 = para + E[,i]

    p1 = obsloglik(paraEi1,x0,delta,t,tevent,num_class,tau,tolEM)
    
    Me[i,] = (p1-p0)/(h) ####  score function, n by 1
    
  }
  
  
  I = Me%*%t(Me)
  
  varmat = sqrt(diag(solve(I)))
  
  varmat
}