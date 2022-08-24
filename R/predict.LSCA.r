predict.LSCA <- function(object,newdata,timepoint=NULL,...){
  
  p <- predict(object$lplr,newdata,type="response")
  
  num_class <- object$num_class
  covx <- object$covx
  # tevent <- object$tevent
  
  for (l in 1:num_class){
    idx <- num_class+(l-1)*length(covx)
    assign(paste0("lp",l), as.matrix(newdata[,covx]) %*% object$zeta[idx:(idx+length(covx)-1)])
  }
  
  for (l in 1:num_class){
    if (l == 1){
      assign(paste0("exp",l), exp(lp1))
    }else{
      assign(paste0("exp",l), exp(lp1+get(paste0("lp",l))+object$zeta[l-1]))
    }
  }
  
  if (is.null(timepoint)){
    chazt <- object$chaz[sapply(newdata$tildet,function(x) which.max(x <= object$tildet))]
  }else{
    chazt <- object$chaz[max(which(object$tildet <= timepoint))]
  }
  
  s_cls <- list()
  for (l in 1:num_class){
    assign(paste0("s",l), exp(-chazt*get(paste0("exp",l))))
    s_cls[[l]] <- get(paste0("s",l))
  }
  
  s <- 0
  
  for (l in 1:num_class){
    s <- s + p[,l]*get(paste0("s",l))
  }
  
  pred <- list()
  pred$p <- p
  pred$s <- s
  pred$s_cls <- s_cls
  
  return(pred)
  
}