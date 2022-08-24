#' lcphR: Latent Class Proportional Hazards Regression
#'
#'
#' @description Conduct latent class analysis with time-to-event data.
#' @param dat Input data.
#' @param num_class Number of latent classes.
#' @param covx Column names of baseline covariates.
#' @param tolEM Stopping criteria for the EM algorithm.
#' @param tolNewton Stopping criteria for the Newton-Raphson algorithm in the M-step.
#' @param maxiterEM Number of maximum iterations of the EM algorithm.
#' @param maxiterNewton Number of maximum iterations of the Newton-Raphson algorithm.
#' @param initial Initialization setting, either 'kmeans', 'null', 'random'. Ignored if probability vectors are provided in 'init.tau'.
#' @param init.tau Provided membership probability based on prior knowledge, a n times c matrix where rows represent subjects and columns represent classes.
#' @param varest TRUE or FALSE, indicating whether conducting variance estimation or not.
#' @param numvar If varest is TRUE, this indicates the numerical method for profile likelihood based variance estimation, either "1st", "1st_rich", or "2nd", indicating 1st order numerical differentiation, 1st order differentiation based on Richardson extrapolation, and 2nd order numerical differentiation, respectively.
#' @param traceplot If TRUE, plot observed data log likelihood and complete data log likelihood versus number of iterations.
#' @param verbose If TRUE, print notifications when point estimation or variance estimation starts.
#' @param warmup Number of warm-up iterations needed for non-informative initial value. Default is 200.
#' @return A list with point estimates (alpha, zeta), variance estimates [ASE (observed log likelihood based), nASESr (profile log likelihood based)].
#' @author Teng Fei. Email: feit1@mskcc.org
#' @references Fei, T., Hanfelt, J., & Peng, L. (2022). Latent Class Analysis with Semi-parametric Proportional Hazards Submodel for Time-to-event Data. arXiv preprint arXiv:2202.00775.
#' @examples
#'
#' set.seed(21227)
#' dat <- simulation.lcphR(n=1000)
#' lcfit <- lcphR(dat,num_class=2,covx=c('Xcov1','Xcov2'),tolEM=1e-3,maxiterEM=100,varest=TRUE,traceplot=TRUE,initial='kmeans')
#'
#' @import VGAM Rcpp RcppArmadillo
#' @useDynLib LCAtoolbox
#' @export

lcphR <- function(dat,num_class,covx,tolEM=1e-7,tolNewton=1e-8,maxiterEM=1000,maxiterNewton=50,initial='null',init.tau=NULL,varest=FALSE,numvar="1st",traceplot=FALSE,verbose=TRUE,warmup=200){

  start = proc.time()[3]
  if (verbose) print('Point estimation started. It may take several minutes till convergence.')

  # require(progress)
  # pb <- progress_bar$new(total=maxiterEM)

  n=length(unique(dat$id))

  dat <- dat[order(dat$tildet,dat$id,dat$time),]

  # extract baseline data
  baseline <- dat[match(unique(dat$id), dat$id),]

  x <- as.matrix(baseline[,covx])
  delta <- baseline$delta
  t <- baseline$tildet
  tevent <- baseline$tildet[as.logical(baseline$delta)]

  lab_class = rep(1,nrow(baseline)*num_class)
  for (i in 2:num_class){
    lab_class[((i-1)*nrow(baseline)+1):((i)*nrow(baseline))] = i
  }

  # Initialize tau0
  if (!is.null(init.tau)){
    tau0 = init.tau
  }else if (initial == 'null'){
    tau0 <- matrix(1/num_class,nrow=n,ncol=num_class)
  }else if (initial == 'kmeans'){
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    # idx <- which(delta==1)
    km = kmeans(t,num_class)$cluster
    for(i in 1:n){
      if(delta[i] == 1){
        tau0[i,km[i]] = 1
      }
    }
    tau0[tau0 < 1e-8] = 1e-8
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
  }else if(initial == 'true'){ # remained for simulation verification
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    for (i in 1:n){
      tau0[i,baseline$latent[i]] = 1
    }
    tau0[tau0<1e-8]=1e-8
  }else if(initial == 'perturbed'){ # remained for simulation verification
    perturbation=0.05
    tau0 <- matrix(0,nrow=n,ncol=num_class)
    u <- runif(n,0,perturbation)
    for (i in 1:n){
      tau0[i,baseline$latent[i]] = 1-u[i]
      tau0[i,-baseline$latent[i]] = u[i]/num_class
    }
  }

  ## EM algorithm

  zeta = rep(0,num_class*ncol(x)+num_class-1)
  alpha = matrix(0,ncol=num_class-1,nrow=ncol(x)+1,byrow=T)
  d = rep(1/length(tevent),length(tevent))
  # taures = postweight(alpha,zeta,x,delta,t,tevent,num_class,d)
  # tau0 = taures$tau

  vars <- paste(covx,collapse="+")
  regression <- paste0("as.factor(class)", " ~ ", vars)

  for (i in 1:maxiterEM){

    alpha0 = alpha
    zeta0 = zeta
    d0 = d
    # pb$tick()
    # Sys.sleep(1/maxiterEM)

    lplr <- suppressWarnings(VGAM::vglm(as.formula(regression),family = multinomial(refLevel = 1),
         weight=tau0,
         data=data.frame(do.call("rbind", rep(list(baseline), num_class)),
                         class = lab_class, tau0 = as.vector(tau0)
         )
         ))

    alpha <- coef(lplr)
    alpha = matrix(alpha,ncol=num_class-1,nrow=ncol(x)+1,byrow=T)

    zeta = newton(maxiterNewton,tolNewton,alpha,zeta,x,delta,t,tevent,num_class,tau0,lambda_EM,loglik_EM)


    lambda_res <- lambda_EM(zeta,x,t,tevent,num_class,tau0)
    taures = postweight(alpha,zeta,x,delta,t,tevent,num_class,lambda_res$d)
    tau1 = taures$tau
    res <- loglik_EM(alpha,zeta,x,delta,t,tevent,num_class,lambda_res$d,lambda_res$d1,lambda_res$d2,tau1)

    chaz = cumsum(lambda_res$d)
    diffPAR = sum(abs(alpha-alpha0),abs(zeta-zeta0),abs(lambda_res$d-d0))

    if(i == 1){
      l0 = res$obsloglik
    }else if(i == 2){
      diffl = res$obsloglik - l0
      l0 = res$obsloglik
    }else if(i == 3){
      a = (res$obsloglik - l0)/(diffl+1e-8)
      diffl = res$obsloglik - l0
      lA = l0 + diffl/(1-a)
      l0 = res$obsloglik
    }else if(i >= 4){
      a = (res$obsloglik - l0)/(diffl+1e-8)
      diffl = res$obsloglik - l0
      difflA = l0 + diffl/(1-a) - lA
      lA = l0 + diffl/(1-a)
      l0 = res$obsloglik

      if ((initial %in% c("true","kmeans","random") & abs(difflA) < tolEM) | i == maxiterEM){
        # Variance estimation
        if (varest){
          if (verbose) print('Numerical variance estimation started.')
          para = c(as.vector(alpha),zeta)
          lambda_res$d = rep(1/length(tevent),length(tevent))
          if (numvar=="1st"){
            nASESr = VarEst_num(para,x0=x,delta=delta,t=t,tevent=tevent,num_class=num_class,tau=tau0,lambda_res,tolEM=0.001)
          }else if(numvar=="1st_rich"){
            nASESr = VarEst_num_richardson(para,x0=x,delta=delta,t=t,tevent=tevent,num_class=num_class,tau=tau0,lambda_res,tolEM=0.001)
          }else if(numvar=="2nd"){
            nASESr = VarEst_num_sec(para,x0=x,delta=delta,t=t,tevent=tevent,num_class=num_class,tau=tau0,lambda_res,tolEM=0.001)
          }
          lambda_res <- lambda_EM(zeta,x,t,tevent,num_class,tau0)
        }
        break
      }else if ((initial=="null" & i > warmup & abs(difflA) < tolEM) | i == maxiterEM){
        # Variance estimation
        if (varest){
          if (verbose) print('Numerical variance estimation started.')
          para = c(as.vector(alpha),zeta)
          lambda_res$d = rep(1/length(tevent),length(tevent))
          if (numvar=="1st"){
            nASESr = VarEst_num(para,x0=x,delta=delta,t=t,tevent=tevent,num_class=num_class,tau=tau0,lambda_res,tolEM=0.001)
          }else if(numvar=="1st_rich"){
            nASESr = VarEst_num_richardson(para,x0=x,delta=delta,t=t,tevent=tevent,num_class=num_class,tau=tau0,lambda_res,tolEM=0.001)
          }else if(numvar=="2nd"){
            nASESr = VarEst_num_sec(para,x0=x,delta=delta,t=t,tevent=tevent,num_class=num_class,tau=tau0,lambda_res,tolEM=0.001)
          }
          lambda_res <- lambda_EM(zeta,x,t,tevent,num_class,tau0)
        }
        break
      }
    }

    d = lambda_res$d
    diffEM = sum((tau1 - tau0)^2)
    tau0 = tau1

    # if (verbose==T){
    #cat(paste('iteration: ',i,'\n','convergence criteria: ',diffPAR,'\n','loglik: ',res$obsloglik,'\n',sep=''))
    # }

    if (traceplot){
      if (is.null(baseline$latent)){
        color <- apply(tau0,1,which.max)
      }else{
        color <- baseline$latent
      }

      pch <- ifelse(delta,4,1)

      par(mfrow=c(1,num_class+1))
      for (l in 1:(num_class)){
        suppressWarnings(plot(t,tau0[1:n,l],ylim=c(0,1),col=color,pch=pch,xlab="event time",ylab=paste0("Posterior Prob of class ",l)))
      }
      suppressWarnings(plot.new())
      text(x = 0.5, y = 0.9, paste("Iteration:",i),
           cex = 2, col = "black")

      text(x = 0.5, y = 0.7, paste("M-step loglik:",round(res$loglik,1)),
           cex = 2, col = "black")

      text(x = 0.5, y = 0.5, paste("Obs Loglik:",round(res$obsloglik,1)),
           cex = 2, col = "black")

      if (i >=4){
        text(x = 0.5, y = 0.3, paste("Convergence:",round(abs(difflA),5)),
           cex = 2, col = "black")
      }

      legend("bottom",c(paste("Class",1:num_class),"Uncensored","Censored"),col=c(1:num_class,1,1),pch=c(rep(1,num_class),4,1))
      # printplot(count=i,res=res,maxiter=maxiterEM)

      # Sys.sleep(1/maxiterEM)

    }

  }
  count = i

  p = taures$p
  res <- loglik_EM(alpha,zeta,x,delta,t,tevent,num_class,lambda_res$d,lambda_res$d1,lambda_res$d2,tau1)
  AIC = -2*res$obsloglik + 2*(length(alpha)+length(zeta))
  BIC = -2*res$obsloglik + log(n)*(length(alpha)+length(zeta))
  ICLBIC = BIC - 2*sum(tau1*log(tau1))

  chaz <- chaz_EM(lambda_res,t,tevent)

  if(varest){
    if (verbose) print('Analytical variance estimation started.')
    target2 <- max(which(tevent < 2))
    target3 <- max(which(tevent < 3))
    target4 <- max(which(tevent < 4))
    ASE = varestcpp(zeta,delta,chaz,tau1,p,x,tevent)
    ASE_chaz = ASE$I
    ASE_par = ASE_chaz[(1:(length(alpha)+length(zeta))),(1:(length(alpha)+length(zeta)))]
    ASE_chaz = ASE_chaz[-(1:(length(alpha)+length(zeta))),-(1:(length(alpha)+length(zeta)))]
    ASEchaz2 = sqrt(rep(1,target2) %*% ASE_chaz[1:target2,1:target2] %*% rep(1,target2))
    ASEchaz3 = sqrt(rep(1,target3) %*% ASE_chaz[1:target3,1:target3] %*% rep(1,target3))
    ASEchaz4 = sqrt(rep(1,target4) %*% ASE_chaz[1:target4,1:target4] %*% rep(1,target4))

    for (c in 1:(num_class-1)){
      assign(paste0("ASEchaz2",c+1), sqrt(c(chaz[target2]*exp(zeta[c]),rep(exp(zeta[c]),target2)) %*% ASE$I[c(length(alpha)+c,(length(alpha)+length(zeta)+1):(length(alpha)+length(zeta)+target2)),c(length(alpha)+c,(length(alpha)+length(zeta)+1):(length(alpha)+length(zeta)+target2))] %*% c(chaz[target2]*exp(zeta[c]),rep(exp(zeta[c]),target2))))
      assign(paste0("ASEchaz3",c+1), sqrt(c(chaz[target3]*exp(zeta[c]),rep(exp(zeta[c]),target3)) %*% ASE$I[c(length(alpha)+c,(length(alpha)+length(zeta)+1):(length(alpha)+length(zeta)+target3)),c(length(alpha)+c,(length(alpha)+length(zeta)+1):(length(alpha)+length(zeta)+target3))] %*% c(chaz[target3]*exp(zeta[c]),rep(exp(zeta[c]),target3))))
      assign(paste0("ASEchaz4",c+1), sqrt(c(chaz[target4]*exp(zeta[c]),rep(exp(zeta[c]),target4)) %*% ASE$I[c(length(alpha)+c,(length(alpha)+length(zeta)+1):(length(alpha)+length(zeta)+target4)),c(length(alpha)+c,(length(alpha)+length(zeta)+1):(length(alpha)+length(zeta)+target4))] %*% c(chaz[target4]*exp(zeta[c]),rep(exp(zeta[c]),target4))))
    }
    ASE$I = NULL

  }

  end = proc.time()[3]
  timediff = end-start
  entropy = 1 + sum(tau1*log(tau1))/(n*log(num_class))
  censor = 1-mean(delta)

  output <- list()
  output$alpha = alpha
  output$zeta = zeta
  output$tildet = baseline$tildet

  if(varest){
    output$ASE = ASE
    # output$Ipar = ASE_par
    # output$nASEI = nASEI
    # output$nASEQ = nASEQ
    # output$nASES = nASES
    output$nASESr = nASESr$ASE
    # output$nASESrIpar = nASESr$I
    target2 <- max(which(t < 2))
    target3 <- max(which(t < 3))
    target4 <- max(which(t < 4))
    output$chaztgt = chaz[c(target2,target3,target4)]
    output$chaztgtASE = c(ASEchaz2,ASEchaz3,ASEchaz4)
    # output$chaztgtASE2 = c(ASEchaz22,ASEchaz32,ASEchaz42)
    # if (num_class > 2){
    #   output$chaztgtASE3 = c(ASEchaz23,ASEchaz33,ASEchaz43)
    # }

  }

  output$tau = tau1
  # output$p = p
  output$loglik = res$loglik
  output$obsloglik = res$obsloglik
  output$AIC = AIC
  output$BIC = BIC
  output$CEBIC = ICLBIC
  output$timediff = timediff
  output$diffEM = diffEM
  output$diffPAR = diffPAR
  output$difflA = difflA
  output$numiter = count
  output$entropy = entropy
  output$censor = censor
  # if (moredetails == T){
    output$chaz = chaz
    output$lplr = lplr
    output$num_class = num_class
    output$covx = covx
  # }

  class(output) <- "LSCA"
  return(output)
}
