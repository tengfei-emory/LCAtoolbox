printplot <- function(count,res,maxiter) {                        
  #print(BIC)                
  if (count == 1){
    plot(count,res$loglik, xlim=c(0,maxiter), type='b', xlab='number of iteration', ylab='loglik')#,ylim=c(res$loglik-1000,res$loglik+1000))
    points(count,res$obsloglik, xlim=c(0,maxiter), type='b',col='red')
  }else if (count > 1){
    points(count,res$loglik, xlim=c(0,maxiter), type='b')  
    points(count,res$obsloglik, xlim=c(0,maxiter), type='b',col='red')
  }
  Sys.sleep(1/maxiter)                   
  #return(BIC)
}