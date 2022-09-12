#' SLTCA: Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood
#'
#' @description Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood
#' @param dat Input data.
#' @param num_class Number of latent classes.
#' @param covx Names of baseline covariates.
#' @param vary Names longitudinal outcomes.
#' @param covgee Names of covariates in GEE models.
#' @param tolEM Stopping criteria for the EM algorithm.
#' @param maxiterEM Number of maximum iterations of the EM algorithm.
#' @param initial Initialization setting, either 'kmcov', 'null', 'random', 'true' or provided in 'init.tau'.
#' @param cor Working correlation structure, currently only supporting "ar1".
#' @param init.tau Provided membership probability based on prior knowledge.
#' @param varest Variance estimation indicator.
#' @param lbound Constant bound of linear projection of approximated likelihood ratio.
#' @param verbose Indicating whether output progress.
#' @param stop.rule Stopping rule, either "PAR" for parameter L2 norm, or "tau" for posterior probability.
#' @return Returns parameter estimation and variance estimation.
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Hart, K. R., Fei, T., & Hanfelt, J. J. (2021). Scalable and robust latent trajectory class analysis using artificial likelihood. Biometrics, 77(3), 1118-1128.
#' @examples 
#' 
#' dat <- simulation.SLTCA(500)
#' res <- SLTCA(dat,num_class=2,covx="baselinecov",vary=paste("y.",1:6,sep=''),
#'              Y_dist=c('poi','poi','bin','bin','normal','normal'),varest=TRUE,verbose=TRUE,stop.rule="tau")
#' 
#' @useDynLib LCAtoolbox, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import VGAM
#' @import geepack
#' @export

SLTCA <- function(dat,num_class,covx,vary,covgee=NULL,Y_dist,tolEM=1e-3,maxiterEM=500,
                  initial='null',cor="ar1",init.tau=NULL,varest=FALSE,lbound=2,verbose=FALSE,
                  stop.rule="PAR",keep.switch=FALSE){
  
  start = proc.time()[3]
  
  if (verbose) {
    print('Point estimation started:')
  }
  
  n=length(unique(dat$id))
  num_feature <- length(Y_dist)
  
  dat <- dat[order(dat$id,dat$time),]
  
  # if (is.null(covgee)){
  #   covgee = covx
  # }
  
  # timedepvars <- paste(c(covx,vary),collapse="+")
  # regression <- paste0("Surv(time,time1,delta)", " ~ ", timedepvars)
  
  # model1 <- coxph(as.formula(regression), data = dat[dat$latent==1,], control = coxph.control(timefix = FALSE))
  # model2 <- coxph(as.formula(regression), data = dat[dat$latent==2,], control = coxph.control(timefix = FALSE))
  # coef1 <- model1$coefficients
  # coef2 <- model2$coefficients
  
  # extract baseline data
  baseline <- dat[match(unique(dat$id), dat$id),]
  baseline$newid <- 1:nrow(baseline)
  
  for (i in 1:nrow(dat)){
    dat$newid[i] = baseline$newid[baseline$id == dat$id[i]]
  }
  
  if (!is.null(covx)){
    x <- as.matrix(baseline[,covx])
  }
  
  lab_class = rep(1,nrow(baseline)*num_class)
  for (i in 2:num_class){
    lab_class[((i-1)*nrow(baseline)+1):((i)*nrow(baseline))] = i
  }
  
  # Initialize tau0
  if (!is.null(init.tau)){
    tau0 = init.tau
  }else if (initial == 'null'){
    tau0 <- matrix(1/num_class,nrow=n,ncol=num_class)
  }else if(initial == 'kmcov'){
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    # idx <- which(delta==1)
    km = kmeans(baseline[,covx],num_class)$cluster
    for(i in 1:n){
      tau0[i,km[i]] = 1
    }
    tau0[tau0 < 1e-8] = 1e-8
  }else if(initial == 'random'){
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    prop <- 0.1
    u <- runif(n,0,1)
    for (i in 1:n){
      if (u[i] < prop){
        tau0[i,ceiling(num_class*runif(1))] = 1
      }
    }
    tau0[tau0 < 1e-8] = 1e-8
  }else if(initial == 'true'){
    tau0 <- matrix(0,nrow=n,ncol=num_class)
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
    tau0[tau0<1e-8]=1e-8
  }else if(initial == 'perturbed'){
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    u <- runif(n,0,0.0001)
    for (i in 1:n){
      tau0[i,baseline$latent[i]] = 1-u[i]
      tau0[i,-baseline$latent[i]] = u[i]/(num_class-1)
    }
    tau0[tau0<1e-8]=1e-8
  }
  
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
  for (i in 1:(length(covgee)+2)){
    beta[[i]] <- matrix(0,ncol=num_class,nrow=num_feature)
  }
  
  for (i in 1:maxiterEM){
    
    alpha0 = alpha
    beta00 = beta0
    beta10 = beta1
    beta_0 = beta
    phi0 = phi
    gamma0 = gamma
    
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
    eqic = 0
    
    mu <- array(0,dim=c(nrow(dat),num_feature,num_class))
    v <- array(0,dim=c(nrow(dat),num_feature,num_class))
    
    longlik <- matrix(0,ncol=num_class,nrow=nrow(dat))
    # plik <- matrix(0,ncol=num_class,nrow=nrow(dat))
    # survlik2 <- matrix(0,ncol=num_class,nrow=nrow(dat))
    # taulik <- matrix(0,ncol=num_class,nrow=nrow(dat))
    
    vars <- paste(covgee,collapse="+")
    
    gee.object <- list()
    
    if (vars != ""){
      regression <- paste0("yy", " ~ time +", vars)
    }else{
      regression <- "yy ~ time"
    }
    
    for (c in 1:num_class){
      
      tau = rep(0,nrow(dat))
      
      # assign posterior membership probability w.r.t class c to the data
      for (ii in 1:n){
        tau[dat$newid==unique(dat$newid)[ii]] = tau0[ii,c]
      }
      
      gee.object[[c]] <- list()
      
      for (j in 1:num_feature){
        
        # yy is the jth longitudinal marker
        yy <- as.numeric(dat[,vary[j]])
        
        # fit corresponding GEE model with weights tau and AR1 correlation structure
        if(Y_dist[j] == 'normal'){
          geefit <- geepack::geeglm(as.formula(regression),family=gaussian,id=newid,waves=num_obs,weights=tau,data=dat,corstr=cor)#,scale.fix=T,scale.value = 1)
        }else if (Y_dist[j] == 'poi'){
          geefit <- geepack::geeglm(as.formula(regression),family=poisson,id=newid,waves=num_obs,weights=tau,data=dat,corstr=cor)#,scale.fix=T,scale.value = 1)
        }else if (Y_dist[j] == 'bin'){
          geefit <- geepack::geeglm(as.formula(regression),family=binomial,id=newid,waves=num_obs,weights=tau,data=dat,corstr=cor)#,scale.fix=T,scale.value = 1)
        }
        
        gee.object[[c]][[j]] <- geefit
        mu[,j,c] = geefit$fitted.values
        v[,j,c] = geefit$family$variance(geefit$fitted.values)
        
        # obtain point estimates
        
        for (ii in 1:(length(covgee)+2)){
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
        
        # use ic function to obtain eqic
        longlik[,c] = longlik[,c] + ic(geefit,Y_dist[j],yy)
        # print(c)
        # print(Y_dist[j])
        # print(summary(ic(geefit,Y_dist[j],yy)))
        # print("weighted:")
        # print(summary(geeweight[,c]*ic(geefit,Y_dist[j],yy)))
      }
    }
    
    # W <- matrix(0,ncol=num_class,nrow=nrow(dat)) # IPW
    # dW <- matrix(0,ncol=num_class,nrow=nrow(dat)) # dIPW/dzeta
    
    y <- as.matrix(dat[,vary])
    
    if (cor=="ar1"){
      ew <- LinProjCpp(mu,v,t(gamma),t(phi),dat$newid,dat$num_obs,dat$time,y)
    }else if (cor=="ind"){
      ew <- LinProjNew(mu,v,gamma,phi,Y_dist,dat,y)
    }
    
    # restrict the upper and lower bound of ew
    ew[ew>10^lbound] = 10^(lbound)
    ew[ew<10^(-lbound)] = 10^(-lbound)
    
    pew <- p*ew
    pewsum = rowSums(pew)
    tau1 = apply(pew,2,function(x) x/pewsum)
    tau1[tau1<1e-8] = 1e-8
    
    # diffPAR = sum(abs(alpha-alpha0),abs(zeta-zeta0),abs(beta0 - beta00),abs(beta1-beta10))
    diffPAR = sum(abs(alpha-alpha0),abs(unlist(beta)-unlist(beta_0)))
    
    if (initial %in% c("null","random","kmcov")){
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
    
    ASE = varestcpp_sltca(dat$time, as.matrix(baseline[,covx]), as.matrix(dat[,vary]), 
                          tau1, p, Y_dist, dat$newid, mu, length(covgee)+2, as.matrix(baseline[,covgee]), phi, gamma)
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
  names(beta) <- c("Intercept","Time",covgee)
  
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
  # output$rangediff = abs(max(dat$tildet[dat$latent==1]) - max(dat$tildet[dat$latent==2]))
  # output$naive = c(coef1,coef2)
  # output$rangesummary = cbind(summary(dat$tildet[dat$latent==1]),summary(dat$tildet[dat$latent==2]))
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