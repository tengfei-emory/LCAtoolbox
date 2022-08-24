VarEst_num_richardson <- function(para,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM){
  
  n <- length(t)
  
  h = 5*n^(-0.5)  #### tuning parameter
  
  p = length(para)
  
  E = h*diag(rep(1,p))  ### p by p
  
  Me = matrix(0,(p),n)
  
  for (i in 1:p){
    
    paraEi1 = para - 2*E[,i]
    paraEi2 = para - E[,i]
    paraEi3 = para + E[,i]
    paraEi4 = para + 2*E[,i]
    
    p1 = obsloglik(paraEi1,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    p2 = obsloglik(paraEi2,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    p3 = obsloglik(paraEi3,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    p4 = obsloglik(paraEi4,x0,delta,t,tevent,num_class,tau,lambda_res,tolEM)
    
    Me[i,] = ( p1 - 8*p2 + 8*p3 - p4 )/(12*h)
    
  }
  
  
  I = Me%*%t(Me)
  
  varmat = sqrt(diag(solve(I)))
  Ipar = solve(I)
  
  output = list(ASE=varmat,I=Ipar)
  return(output)
}