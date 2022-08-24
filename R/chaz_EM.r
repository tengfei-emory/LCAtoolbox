chaz_EM <- function(lambda_res,t,tevent){
  
  n <- length(t)
  chaz <- rep(0,n)
  d <- lambda_res$d
  cumd <- cumsum(d)
  
  for (i in 1:length(tevent)){
    if (i < length(tevent)){
      idx <- which(t >= tevent[i] & t < tevent[i+1])
    }else{
      idx <- which(t >= tevent[i])
    }
    chaz[idx] = cumd[i]
  }
  
  chaz
}
