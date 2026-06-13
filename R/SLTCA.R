#' SLTCA: Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood
#'
#' @description Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood
#' @param dat Input data in long format (one row per subject-visit), with columns id, time,
#'   num_obs (visit/wave index), the longitudinal outcomes and covariates. Outcomes may be
#'   missing (NA) at any subset of visits: each outcome enters the model through the visits
#'   at which it is observed, by conditional independence of the outcomes given the class.
#' @param num_class Number of latent classes.
#' @param covx Names of baseline covariates.
#' @param vary Names of longitudinal outcomes.
#' @param covgee Names of covariates in GEE models.
#' @param Y_dist Distributions of longitudinal outcomes, a vector with elements 'normal', 'bin' or 'poi'.
#' @param tolEM Stopping criteria for the EM algorithm.
#' @param maxiterEM Number of maximum iterations of the EM algorithm.
#' @param initial Initialization setting, either 'kmcov', 'null', 'random', 'true' or provided in 'init.tau'.
#' @param cor Working correlation structure: "ind" (independence), "ar1" or "exch" (exchangeable).
#' @param ipw A vector of inverse probability weight to account for informative censoring.
#' @param init.tau Provided membership probability based on prior knowledge.
#' @param varest Variance estimation indicator.
#' @param lbound Constant bound of linear projection of approximated likelihood ratio.
#' @param verbose Indicating whether output progress.
#' @param stop.rule Stopping rule, either "PAR" for parameter L2 norm, or "tau" for posterior probability.
#' @param keep.switch Whether subjects with near-degenerate membership keep their class between iterations.
#' @return Returns parameter estimation and variance estimation.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Hart, K. R., Fei, T., & Hanfelt, J. J. (2021). Scalable and robust latent trajectory class analysis using artificial likelihood. Biometrics, 77(3), 1118-1128.
#' @examples
#'
#' dat <- simulation.SLTCA(500)
#' res <- SLTCA(dat,num_class=2,covx="baselinecov",vary=paste("y.",1:6,sep=''),covgee="time",
#'              Y_dist=c('poi','poi','bin','bin','normal','normal'),varest=TRUE,verbose=TRUE,stop.rule="tau")
#'
#' @useDynLib LCAtoolbox, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import VGAM
#' @import geepack
#' @export

SLTCA <- function(dat,num_class,covx,vary,covgee,Y_dist,tolEM=1e-3,maxiterEM=500,
                  initial='null',cor="ar1",ipw=1,init.tau=NULL,varest=FALSE,lbound=2,verbose=FALSE,
                  stop.rule="PAR",keep.switch=FALSE){

  start = proc.time()[3]

  # map working correlation to geepack's corstr; the same short code is used in the C++ routines
  corstr <- switch(cor,
                   ind = "independence",
                   ar1 = "ar1",
                   exch = "exchangeable",
                   stop("'cor' must be one of \"ind\", \"ar1\" or \"exch\""))

  if (verbose) {
    print('Point estimation started:')
  }

  n = length(unique(dat$id))
  num_feature <- length(Y_dist)

  dat <- dat[order(dat$id,dat$time),]
  dat$newid <- match(dat$id, unique(dat$id))

  # extract baseline data
  baseline <- dat[match(unique(dat$id), dat$id),]
  baseline$newid <- 1:nrow(baseline)

  if (!is.null(covx)){
    x <- as.matrix(baseline[,covx])
  }

  # observation indicators: which visits have each outcome observed
  yobs <- !is.na(as.matrix(dat[,vary,drop=FALSE]))
  if (any(colSums(yobs) == 0)){
    stop("Outcome(s) ", paste(vary[colSums(yobs)==0],collapse=", "), " have no observed values")
  }

  if (length(ipw) == 1) ipw <- rep(ipw,nrow(dat))

  lab_class = rep(1:num_class,each=nrow(baseline))

  tau0 <- init_tau(init.tau,initial,num_class,n,baseline,covx)

  ## EM algorithm

  if (is.null(covx)){
    alpha = matrix(0,ncol=num_class-1,nrow=1,byrow=T)
  }else{
    alpha = matrix(0,ncol=num_class-1,nrow=ncol(x)+1,byrow=T)
  }

  beta0 <- matrix(0,ncol=num_class,nrow=num_feature)
  beta1 <- matrix(0,ncol=num_class,nrow=num_feature)
  phi <- matrix(0,ncol=num_class,nrow=num_feature)
  gamma <- matrix(0,ncol=num_class,nrow=num_feature)
  beta = list()
  for (i in 1:(length(covgee)+1)){
    beta[[i]] <- matrix(0,ncol=num_class,nrow=num_feature)
  }

  geefamily <- list(normal = gaussian, poi = poisson, bin = binomial)

  for (i in 1:maxiterEM){

    alpha0 = alpha
    beta_0 = beta

    if (!is.null(covx)){
      vars <- paste(covx,collapse="+")
      regression <- paste0("as.factor(class)", " ~ ", vars)
    }else{
      regression <- "as.factor(class) ~ 1"
    }

    if (num_class > 1){
      lcfit <- vglm(as.formula(regression),family = multinomial(refLevel = 1),
                    weight=tau0,
                    data=data.frame(do.call("rbind", rep(list(baseline), num_class)),
                                    class = lab_class, tau0 = as.vector(tau0)
                    ))
      alpha <- coef(lcfit)
      alpha = matrix(alpha,ncol=ncol(alpha0),nrow=nrow(alpha0),byrow=T)
      p <- as.matrix(lcfit@fitted.values[1:n,])
    }else if (num_class==1){
      p <- as.matrix(rep(1,n),nrow=n,ncol=1)
    }

    # GEE

    mu <- array(0,dim=c(nrow(dat),num_feature,num_class))
    v <- array(0,dim=c(nrow(dat),num_feature,num_class))

    longlik <- matrix(0,ncol=num_class,nrow=nrow(dat))

    vars <- paste(covgee,collapse="+")

    gee.object <- list()

    if (vars != ""){
      regression <- paste0("yy", " ~ ", vars)
    }else{
      regression <- "yy ~ time"
    }

    for (c in 1:num_class){

      # assign posterior membership probability w.r.t class c to the data
      tau <- tau0[dat$newid,c]

      gee.object[[c]] <- list()

      for (j in 1:num_feature){

        # fit the GEE model for the jth outcome on the visits where it is observed,
        # with weights tau and the requested working correlation structure
        obsj <- yobs[,j]
        datj <- dat[obsj,]
        datj$yy <- as.numeric(datj[,vary[j]])
        datj$geeweight <- (tau/ipw)[obsj]

        geefit <- geepack::geeglm(as.formula(regression),family=geefamily[[Y_dist[j]]],
                                  id=newid,waves=num_obs,weights=geeweight,
                                  data=datj,corstr=corstr)

        gee.object[[c]][[j]] <- geefit

        # fitted means and variance functions at all visits (used in the linear
        # projection; visits with missing y are skipped there via NA in y)
        mu[,j,c] = predict(geefit,newdata=dat,type="response")
        v[,j,c] = geefit$family$variance(mu[,j,c])

        # obtain point estimates
        for (ii in 1:(length(covgee)+1)){
          beta[[ii]][j,c] <- coef(geefit)[ii]
        }

        beta0[j,c] <- coef(geefit)[1]
        beta1[j,c] <- coef(geefit)[2]

        phi[j,c] <- geefit$geese$gamma
        if (cor == "ind"){
          gamma[j,c] = 1
        }else{
          gamma[j,c] <- geefit$geese$alpha
        }

        # use ic function to obtain eqic; observed visits only
        longlik[obsj,c] = longlik[obsj,c] + ic(geefit,Y_dist[j],datj$yy)
      }
    }

    y <- as.matrix(dat[,vary])

    ew <- LinProjCpp(mu,v,t(gamma),t(phi),dat$newid,dat$num_obs,dat$time,y,cor)

    # restrict the upper and lower bound of ew
    ew[ew>10^lbound] = 10^(lbound)
    ew[ew<10^(-lbound)] = 10^(-lbound)

    pew <- p*ew
    pewsum = rowSums(pew)
    tau1 = apply(pew,2,function(x) x/pewsum)
    tau1[tau1<1e-8] = 1e-8

    diffPAR = sum(abs(alpha-alpha0),abs(unlist(beta)-unlist(beta_0)))

    if (initial %in% c("null","random","kmcov")){
      if (keep.switch){
        tau01 = tau0
        for (idx in 1:n){
          if (sum(tau01[idx,] > 0.995) == 0){
            tau01[idx,] <- tau1[idx,]
          }
        }
        rownames(tau01) = rownames(tau1)
      }else{
        tau01 = tau1
      }
      diffEM = sum((tau01 - tau0)^2)/num_class
      tau0 = tau01
      tau1 = tau01
    }else if (initial=="true"){
      tau01 = tau1
      if (i > 1){
        if (keep.switch){
          tau01 = tau0
          for (idx in 1:n){
            if (sum(tau01[idx,] > 0.995) == 0){
              tau01[idx,] <- tau1[idx,]
            }
          }
        }else{
          tau01 = tau1
        }
      }
      diffEM = sum((tau01 - tau0)^2)/num_class
      tau0 = tau01
      tau1 = tau01
    }

    if (verbose){
      cat(paste('Convergence criteria: ',round(diffEM,5),' ', round(diffPAR,5),'\n',sep=''))
      plot(tau0[,1])
    }

    if (num_class == 1){
      break
    }else if (initial %in% c("null","random","kmcov") | (initial=="true" & num_class < max(baseline$latent))){
      if (max(tau1) > 0.7){
        if (stop.rule == "PAR"){
          if (abs(diffPAR) < tolEM | i == maxiterEM){
            break
          }
        }else if (stop.rule == "tau"){
          if (abs(diffEM) < tolEM | i == maxiterEM){
            break
          }
        }
      }
    }else if (initial=="true"){
      if (stop.rule == "PAR"){
        if (abs(diffPAR) < tolEM | i == maxiterEM){
          break
        }
      }else if (stop.rule == "tau"){
        if (abs(diffEM) < tolEM | i == maxiterEM){
          break
        }
      }
    }

  }

  count=i

  end = proc.time()[3]
  timediff = end-start
  entropy = 1 + sum(tau1*log(tau1))/(n*log(num_class))

  EQIC.long=0

  for (ii in 1:n){
    qloglik <- colSums(matrix(longlik[dat$newid==ii,],ncol=num_class))
    EQIC.long = EQIC.long - 2*sum(qloglik)
  }

  QEAIC.long = EQIC.long + 2*(length(alpha) + length(unlist(beta)))
  QEBIC.long = EQIC.long + log(n)*(length(alpha) + length(unlist(beta)))
  CEEQIC.long = QEBIC.long - 2*sum(tau1*log(tau1))

  modal <- apply(tau1,1,which.max)

  if (varest == TRUE){

    # rows with missing y (NA -> NaN on the C++ side) are skipped per outcome
    if (is.null(covx)){
      ASE = varestcpp_sltca_prob(dat$time, as.matrix(dat[,vary]),
                                 tau1, p, Y_dist, dat$newid, dat$num_obs, mu, length(covgee)+1, as.matrix(baseline[,covgee]), phi, gamma, cor)
    }else{
      ASE = varestcpp_sltca(dat$time, as.matrix(baseline[,covx]), as.matrix(dat[,vary]),
                            tau1, p, Y_dist, dat$newid, dat$num_obs, mu, length(covgee)+1, as.matrix(baseline[,covgee]), phi, gamma, cor)
    }

    rownames(ASE$betarb) <- rownames(ASE$betai) <- paste(rep(paste("Class",1:num_class),each=length(vary)*(1+length(covgee))),
                                                         rep(vary,num_class*(1+length(covgee))),
                                                         rep(rep(c("Intercept",covgee),each=length(vary)),num_class)
    )

    if (is.null(covx)){
      rownames(ASE$alpharb) <- rownames(ASE$alphai) <- paste("Class",2:num_class)
    }else{
      rownames(ASE$alpharb) <- rownames(ASE$alphai) <- paste(rep(paste("Class",2:num_class),each=length(covx)+1),
                                                             rep(c("Intercept",covx),num_class-1)
      )
    }

  }

  rownames(beta0) <- rownames(beta1) <- rownames(phi) <- rownames(gamma) <- vary
  if (is.null(covx)){
    rownames(alpha) <- "Intercept"
  }else{
    rownames(alpha) <- c("Intercept",covx)
  }

  colnames(beta0) <- colnames(beta1) <- colnames(phi) <- colnames(gamma) <- paste("Class",1:num_class)
  colnames(alpha) <- paste("Class",2:num_class)
  beta <- lapply(beta,function(x) {colnames(x) <- paste("Class",1:num_class);x})
  beta <- lapply(beta,function(x) {rownames(x) <- vary;x})
  names(beta) <- c("Intercept",covgee)

  output <- list()
  output$num_class = num_class
  output$alpha = alpha
  output$beta0 = beta0
  output$beta1 = beta1
  output$beta = beta
  output$phi = phi
  output$gamma = gamma
  output$entropy = entropy
  output$diffEM = diffEM
  output$diffPAR = diffPAR
  output$numiter = count
  output$timediff = timediff
  output$cor = cor
  if (varest == T){
    output$ASE = ASE
  }

  if (!is.null(dat$latent)){
    ARI <- mclust::adjustedRandIndex(baseline$latent,modal)
    output$ARI = ARI
  }

  output$QIC = list(QEAIC.long=QEAIC.long,QEBIC.long=QEBIC.long,CEEQIC.long=CEEQIC.long)
  output$vary = vary
  output$Y_dist = Y_dist
  output$modal = modal

  output$covx = covx
  output$covgee = covgee
  output$tau = tau1

  class(output) <- "SLTCA"

  return(output)
}

# Initialize the posterior membership probability matrix tau0 (n x num_class)
init_tau <- function(init.tau,initial,num_class,n,baseline,covx){

  if (!is.null(init.tau)){
    return(init.tau)
  }

  if (initial == 'null'){
    tau0 <- matrix(1/num_class,nrow=n,ncol=num_class)
    return(tau0)
  }

  tau0 <- matrix(0,nrow=n,ncol=num_class)

  if (initial == 'kmcov'){
    km = kmeans(baseline[,covx],num_class)$cluster
    for(i in 1:n){
      tau0[i,km[i]] = 1
    }
  }else if (initial == 'random'){
    prop <- 0.1
    u <- runif(n,0,1)
    for (i in 1:n){
      if (u[i] < prop){
        tau0[i,ceiling(num_class*runif(1))] = 1
      }
    }
  }else if (initial == 'true'){
    if (max(baseline$latent) == num_class){
      for (i in 1:n){
        tau0[i,baseline$latent[i]] = 1
      }
    }else if (max(baseline$latent) < num_class){
      for (i in 1:n){
        if (baseline$latent[i] < max(baseline$latent)){
          tau0[i,baseline$latent[i]] = 1
        }else{
          tau0[i,sample((max(baseline$latent)):num_class,1)] = 1
        }
      }
    }else{
      for (i in 1:n){
        if (baseline$latent[i] < num_class){
          tau0[i,baseline$latent[i]] = 1
        }else{
          tau0[i,num_class] = 1
        }
      }
    }
  }else if (initial == 'perturbed'){
    u <- runif(n,0,0.0001)
    for (i in 1:n){
      tau0[i,baseline$latent[i]] = 1-u[i]
      tau0[i,-baseline$latent[i]] = u[i]/(num_class-1)
    }
  }else{
    stop("Unknown 'initial' setting: ", initial)
  }

  tau0[tau0 < 1e-8] = 1e-8
  tau0
}
