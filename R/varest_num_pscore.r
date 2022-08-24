VarEst_num_pscore <- function(para,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM){
  
  n <- length(t)
  
  p = length(para)
  
  I <- matrix(0, p, p)
  
  h = 10^(-5)  #### tuning parameter
  
  E = h * diag(rep(1,p))  ### p by p
  
  for (i in 1:p){
    
    muEi1 = para - 2*E[,i]
    muEi2 = para - E[,i]
    muEi3 = para + E[,i]
    muEi4 = para + 2*E[,i]
    
    I1 = pscore(muEi1,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    I2 = pscore(muEi2,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    I3 = pscore(muEi3,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    I4 = pscore(muEi4,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    
    I[i,] = -( I1 - 8*I2 + 8*I3 - I4 )/(12*h)
    
  }
  
  
  I = (I+t(I))/2
  
  sd.theta = sqrt(diag(solve(I)))
  
  sd.theta
  
}