#' predict.SLTCA: Predict membership probabilities from a fitted SLTCA model
#'
#' @description Computes posterior latent class membership probabilities for
#'   (new) subjects using their longitudinal observations up to a target time.
#'   newdata must be in the same long format as the data used to fit the model
#'   (columns id, time, num_obs, the outcomes and covariates). Outcomes may be
#'   missing (NA) at any subset of visits.
#' @param object Fitted SLTCA object.
#' @param newdata New data used for prediction.
#' @param target.time Time at which the prediction should be performed; only
#'   observations with time <= target.time are used. Defaults to Inf (use all).
#' @param ... Not used.
#' @return A list with tau (matrix of posterior membership probabilities, one
#'   row per subject), modal (modal class assignment) and data (the subset of
#'   newdata used for prediction).
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Hart, K. R., Fei, T., & Hanfelt, J. J. (2021). Scalable and robust latent trajectory class analysis using artificial likelihood. Biometrics, 77(3), 1118-1128.
#' @examples
#'
#' dat <- simulation.SLTCA(500)
#' res <- SLTCA(dat,num_class=2,covx="baselinecov",vary=paste("y.",1:6,sep=''),covgee="time",
#'              Y_dist=c('poi','poi','bin','bin','normal','normal'),varest=TRUE,verbose=TRUE,stop.rule="tau")
#' # posterior membership using all observations
#' pred <- predict(res,newdata=dat)
#'
#' # early prediction using only observations up to time 1
#' pred1 <- predict(res,newdata=dat,target.time=1)
#' head(pred1$tau)
#' table(pred1$modal)
#'
#' @method predict SLTCA
#' @export

predict.SLTCA <- function(object,newdata,target.time=Inf,...){

  dat <- newdata[newdata$time <= target.time,,drop=FALSE]
  if (nrow(dat) == 0){
    stop("No observations in 'newdata' with time <= ", target.time)
  }

  dat <- dat[order(dat$id,dat$time),]
  dat$newid <- match(dat$id, unique(dat$id))
  baseline <- dat[match(unique(dat$id), dat$id),]

  num_class <- object$num_class
  vary <- object$vary
  num_feature <- length(vary)

  # prior membership probabilities from the baseline covariate model
  if (is.null(object$covx)){
    xb <- matrix(1,nrow=nrow(baseline),ncol=1)
  }else{
    xb <- as.matrix(cbind(1,baseline[,object$covx]))
  }
  p <- cbind(1,exp(xb %*% object$alpha))
  p <- p/rowSums(p)

  # GEE design matrix: intercept + covgee (SLTCA falls back to ~ time)
  covgee <- object$covgee
  if (is.null(covgee) || length(covgee) == 0) covgee <- "time"
  xmat <- as.matrix(cbind(1,dat[,covgee]))

  linkinv <- list(normal = identity,
                  poi = exp,
                  bin = function(eta) 1/(1+exp(-eta)))
  varfun <- list(normal = function(mu) rep(1,length(mu)),
                 poi = identity,
                 bin = function(mu) mu*(1-mu))

  # fitted means and variance functions from the estimated GEE coefficients
  mu <- array(0,dim=c(nrow(dat),num_feature,num_class))
  v <- array(0,dim=c(nrow(dat),num_feature,num_class))

  for (c in 1:num_class){
    for (j in 1:num_feature){
      betajc <- sapply(object$beta, function(b) b[j,c])
      eta <- xmat %*% betajc
      mu[,j,c] = linkinv[[object$Y_dist[j]]](eta)
      v[,j,c] = varfun[[object$Y_dist[j]]](mu[,j,c])
    }
  }

  y <- as.matrix(dat[,vary])

  ew <- LinProjCpp(mu,v,t(object$gamma),t(object$phi),dat$newid,dat$num_obs,dat$time,y,object$cor)

  lbound <- if (is.null(object$lbound)) 2 else object$lbound
  ew[ew>10^lbound] = 10^(lbound)
  ew[ew<10^(-lbound)] = 10^(-lbound)

  tau <- p*ew/rowSums(p*ew)
  colnames(tau) <- paste0("tau",1:num_class)
  rownames(tau) <- baseline$id

  return(list(tau = tau,
              modal = apply(tau,1,which.max),
              data = dat))

}
