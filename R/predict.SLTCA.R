#' predict.SLTCA: Predict membership probabilities from a fitted SLTCA model 
#'
#' @param fit Fitted SLTCA object.
#' @param newdata New data used for prediction.
#' @param target.time Time at which the prediction should be performed.
#' @return Returns a matrix of latent class membership probabilities.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Hart, K. R., Fei, T., & Hanfelt, J. J. (2021). Scalable and robust latent trajectory class analysis using artificial likelihood. Biometrics, 77(3), 1118-1128.
#' @examples 
#' 
#' dat <- simulation.SLTCA(500)
#' res <- SLTCA(dat,num_class=2,covx="baselinecov",vary=paste("y.",1:6,sep=''),covgee="time",
#'              Y_dist=c('poi','poi','bin','bin','normal','normal'),varest=TRUE,verbose=TRUE,stop.rule="tau")
#' pred <- predict.SLTCA(res,newdata=dat,target.time=0)
#' 
#' @useDynLib LCAtoolbox, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import VGAM
#' @import geepack
#' @export

predict.SLTCA <- function(fit,
                          newdata,
                          target.time){
  
  data_filtered <- newdata[newdata$time <= target.time,]

  n=length(unique(data_filtered$id))
  data_filtered <- data_filtered[order(data_filtered$id,data_filtered$time),]
  rownames(data_filtered) <- 1:nrow(data_filtered)
  
  # extract baseline data_filtereda
  baseline <- data_filtered[match(unique(data_filtered$id), data_filtered$id),]
  baseline$newid <- 1:nrow(baseline)
  
  for (i in 1:nrow(data_filtered)){
    data_filtered$newid[i] = baseline$newid[baseline$id == data_filtered$id[i]]
  }
  
  num_class <- fit$num_class
  num_feature <- nrow(fit$phi)
  vary <- fit$vary
  gee.object <- fit$gee.object
  xbaselinemat <- as.matrix(cbind(1,baseline[,fit$covx]))
  xmat <- as.matrix(cbind(1,data_filtered[,fit$covgee]))
  
  p <- cbind(1,exp(xbaselinemat %*% fit$alpha))
  p <- p/rowSums(p)
  
  mu <- array(0,dim=c(nrow(data_filtered),num_feature,num_class))
  v <- array(0,dim=c(nrow(data_filtered),num_feature,num_class))
  
  for (c in 1:num_class){
    
    for (j in 1:num_feature){
      
      # yy is the jth longitudinal marker
      yy <- as.numeric(data_filtered[,vary[j]])
      
      mu[,j,c] = gee.object[[c]][[j]]$family$linkinv(xmat %*% gee.object[[c]][[j]]$coefficients)
      
      v[,j,c] = gee.object[[c]][[j]]$family$variance(mu[,j,c])
      
    }
  }
  
  gamma <- fit$gamma
  phi <- fit$phi
  cor <- fit$cor
  
  y <- as.matrix(data_filtered[,vary])
  
  ew <- LinProjCpp(mu,v,t(gamma),t(phi),data_filtered$newid,data_filtered$num_obs,data_filtered$time,y,cor)
  
  ew[ew>10^2] = 10^(2)
  ew[ew<10^(-2)] = 10^(-2)
  
  tau <- p*ew/rowSums(p*ew)
  colnames(tau) <- paste0("tau",1:num_class)

  return(list(tau = tau,
              data = data_filtered))
  
}