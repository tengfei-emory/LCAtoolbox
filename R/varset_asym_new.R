VarEst_asym_new <- function(zeta,phi,gamma,Mu,t0,t1,tevent,delta,x,xy,y,bsp,Y_dist,ew,d,tau,p,id,last){
  require(Matrix)
  num_class = ncol(tau)
  num_feature = ncol(y)
  n <- nrow(x)
  
  # zetaxy <- zeta[1:(ncol(xy)*num_class)]
  # zetabsp <- zeta[(ncol(xy)*num_class+1):length(zeta)]
  
  # a: multinomial logistic
  # g: survival coefficient
  # l: baseline hazard jumps
  # t: GEE
  BBa <- 0
  BBag <- 0
  BBal <- 0
  BBat <- 0
  BBg <- 0
  BBgl <- 0
  BBlg <- 0
  BBtg <- 0
  BBl <- 0
  BBtl <- 0
  BBt <- 0
  BBgt <- 0
  BBlt <- 0
  BBta <- 0
  S <- 0
  
  Ba <- 0
  Bag <- 0
  Bal <- 0
  Bat <- list()
  
  Bg <- 0
  Bgl <- Matrix(0,nrow=length(zeta),ncol=length(tevent),sparse=TRUE)
  Bgt <- list()
  
  Bl <- Matrix(0,nrow=length(tevent),ncol=length(tevent),sparse=TRUE)
  Blg <- Matrix(0,nrow=length(tevent),ncol=length(zeta),sparse=TRUE)
  Blt <- list()
  
  Btl <- list()
  Bt <- list()
  Bta <- list()
  Btg <- list()
  
  mu <- list()
  v <- list()
  dmu <- list()
  dv <- list()
  
  for (k in 1:length(tevent)){
    if (k == 1){
      idx <- which(t1 <= tevent[k] & last == 1)
    }else if (k < length(tevent)){
      idx <- which(t1 < tevent[k+1] & t1 >= tevent[k] & last == 1)
    }else{
      idx <- which(t1 >= tevent[k] & last == 1)
    }
    
    ids <- unique(id[idx])
    print(ids)
    
    for (i in ids){
      
      bigBa <- 0
      bigBg <- 0
      bigBl <- 0
      bigBt <- list()
      bigBht <- list()
      
      xi = c(1,x[i,])
      xyi = xy[id==i,]
      t0i = t0[id==i]
      t1i = t1[id==i]
      deltai = delta[id==i]
      # d_diff = sapply(t1[id==i],function(x) sum(d[tevent <= x]))
      # d_diff <- c(d_diff[1],diff(d_diff))
      
      for (c in 1:num_class){
        
        Bt[[c]] <- 0
        Bat[[c]] <- 0
        Bta[[c]] <- 0
        Btg[[c]] <- 0
        Bgt[[c]] <- 0
        Btl[[c]] <- 0
        Blt[[c]] <- 0
        
        # multinomial logistic
        Dalpha <- -outer(p[i,-1],p[i,-1])
        diag(Dalpha) <- p[i,-1]*(1-p[i,-1])
        qalpha <- rep(0,(num_class-1)*(ncol(x)+1))
        for (j in 2:num_class){
          qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] = qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] - p[i,j]*xi
        }
        if (c > 1){
          qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] = qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] + xi
        }
        Ba = Ba + tau[i,c]*kronecker(Dalpha,outer(xi,xi)) - tau[i,c]*qalpha%o%qalpha
        bigBa = bigBa + tau[i,c]*qalpha
        
        # survival
        
        zt <- matrix(0,nrow=sum(id==i),ncol=ncol(xy)*num_class)
        zt[,((c-1)*ncol(xy)+1):(c*ncol(xy))] <- xyi
        bspi <- matrix(0,nrow=k,ncol=ncol(bsp)*(num_class-1))
        if (c > 1){
          bspi[,((c-2)*ncol(bsp)+1):((c-1)*ncol(bsp))] <- bsp[1:k,]
        }
      
        if (max(t1i) >= tevent[k]){
          
          qzeta <- rep(0,length(zeta))
          qlambda <- rep(0,length(tevent))
          
          for (kk in 1:k){
          
            zkk <- c(zt[which(t0i < tevent[kk] & t1i >= tevent[kk]),],bspi[kk,])
            
            qzeta = qzeta - zkk*as.vector(exp(zkk %*% zeta))*d[kk]
            qlambda[kk] = qlambda[kk] - as.vector(exp(zkk %*% zeta))*d[kk]
            
            Bg = Bg + tau[i,c]*outer(zkk,zkk)*as.vector(exp(zkk %*% zeta))*d[kk]
            Bl[kk,kk] = Bl[kk,kk] + tau[i,c]*as.vector(exp(zkk %*% zeta))
            Bgl[,kk] = Bgl[,kk] + tau[i,c]*zkk*as.vector(exp(zkk %*% zeta))
            Blg[kk,] = Blg[kk,] + tau[i,c]*zkk*as.vector(exp(zkk %*% zeta))*d[kk]

          }
          
          qzeta <- qzeta + zkk*deltai[length(deltai)]
          qlambda[k] <- qlambda[k] + deltai[length(deltai)]
          
          Bg = Bg - tau[i,c]*qzeta%o%qzeta
          Bl = Bl - tau[i,c]*qlambda%o%qlambda
          Bgl = Bgl - tau[i,c]*qzeta%o%qlambda
          Blg = Blg - tau[i,c]*qlambda%o%qzeta
          
          bigBg = tau[i,c]*qzeta
          bigBl = tau[i,c]*qlambda
        }
        
        
        # GEE
        
        # longitudinal part
        
        mu[[c]] <- list()
        v[[c]] <- list()
        dmu[[c]] <-list()
        dv[[c]] <- list()
        Wi = as.vector(W[dat$newid==i,c])
        dWi = as.vector(dW[dat$newid==i,c])
        
        for (j in 1:num_feature){
          if (Y_dist[j] == 'normal'){
            mulink <- function(x) x
            varlink <- function(x) rep(1,length(x))
            dmulink <- function(x) 1
            dvarlink <- function(x) 0
          }else if (Y_dist[j] == 'poi'){
            mulink <- function(x) exp(x)
            varlink <- function(x) x
            dmulink <- function(x) exp(x)
            dvarlink <- function(x) 1
          }else if (Y_dist[j] == 'bin'){
            mulink <- function(x) exp(x)/(1+exp(x))
            varlink <- function(x) x*(1-x)
            dmulink <- function(x) exp(x)/(1+exp(x))^2 #(exp(x)*(1+exp(x)) - exp(x)*exp(x))/(1+exp(x))^2
            dvarlink <- function(x) 1-2*x
          }
          
          mu[[c]][[j]] <- Mu[dat$newid==i,j,c]
          if (length(mu[[c]][[j]]) == 1){
            v[[c]][[j]] <- phi[j,c]*diag(sqrt(varlink(mu[[c]][[j]])),nrow=1) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(sqrt(varlink(mu[[c]][[j]])),nrow=1)
          }else{
            v[[c]][[j]] <- phi[j,c]*diag(as.vector(sqrt(varlink(mu[[c]][[j]])))) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j]]))))
          }
          dmu[[c]][[j]] <- cbind(dmulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$newid==i]))),
                                 dmulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$newid==i])))*dat$time[dat$newid==i])
          
          if (length(mu[[c]][[j]]) == 1){
            dv_uv <- phi[j,c]*diag(dvarlink(mu[[c]][[j]])/(2*sqrt(varlink(mu[[c]][[j]]))),nrow=1) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(sqrt(varlink(mu[[c]][[j]])),nrow=1)
            dv_uv_dmu <- outer(dv_uv,dmu[[c]][[j]])
            dv[[c]][[j]] <- array(0,dim=c(length(mu[[c]][[j]]),length(mu[[c]][[j]]),2)) # here I only consider intercept and time as the covariates (2 dimensional)
            dv[[c]][[j]][1,1,] = (dv_uv_dmu[1,1,1,] + dv_uv_dmu[1,1,1,])
          }else{
            dv_uv <- phi[j,c]*diag(as.vector(dvarlink(mu[[c]][[j]])/(2*sqrt(varlink(mu[[c]][[j]]))))) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j]]))))
            dv_uv_dmu <- outer(dv_uv,dmu[[c]][[j]])
            dv[[c]][[j]] <- array(0,dim=c(length(mu[[c]][[j]]),length(mu[[c]][[j]]),2)) # here I only consider intercept and time as the covariates (2 dimensional)
            for (uu in 1:length(mu[[c]][[j]])){
              for (vv in 1:length(mu[[c]][[j]])){
                dv[[c]][[j]][uu,vv,] = (dv_uv_dmu[uu,vv,uu,] + dv_uv_dmu[vv,uu,vv,])
              }
            }
          }
        }
        
        muc <- unlist(mu[[c]])
        mu01 <- unlist(mu[[1]])
        dmuc <- dmu[[c]]
        vc <- v[[c]]
        dmu1 <- dmu[[1]]
        v1 <- v[[1]]
        
        dv1_int <- Matrix::bdiag(lapply(dv[[1]],function(x) (x[,,1])))
        dv1_t <- Matrix::bdiag(lapply(dv[[1]],function(x) (x[,,2])))
        dvc_int <- Matrix::bdiag(lapply(dv[[c]],function(x) (x[,,1])))
        dvc_t <- Matrix::bdiag(lapply(dv[[c]],function(x) (x[,,2])))
        vinvc <- as.matrix(solve(Matrix::bdiag(vc)))
        vinv1 <- as.matrix(solve(Matrix::bdiag(v1)))
        
        qt <- as.vector(t(as.matrix(Matrix::bdiag(dmuc))) %*% vinvc %*% ((as.vector(yi)-muc)/as.vector(Wi)) )
        
        if (c == 1){
          ht = rep(0,num_feature*2)
        }else{
          ht1 <- as.vector(t(as.matrix(Matrix::bdiag(dmuc)) - as.matrix(Matrix::bdiag(dmu1))) %*% (vinvc %*% ((as.vector(yi)-muc)) + vinv1 %*% ((as.vector(yi)-mu01))))/(2*log(ew[i,c]))
          ht2_int <- (mu01 - muc) * (- vinvc %*% dvc_int %*% vinvc %*% ((as.vector(yi)-muc))
                                     - vinv1 %*% dv1_int %*% vinv1 %*% ((as.vector(yi)-mu01)))/(2*log(ew[i,c]))
          ht2_t <- (mu01 - muc) * (- vinvc %*% dvc_t %*% vinvc %*% ((as.vector(yi)-muc))
                                   - vinv1 %*% dv1_t %*% vinv1 %*% ((as.vector(yi)-mu01)))/(2*log(ew[i,c]))
          transmat <- kronecker(diag(6),rep(1,length(mu[[c]][[j]])))
          # print(i)
          # print(ht2_int)
          # print(ht2_t)
          ht2 = as.vector(t(cbind(as.matrix(ht2_int),as.matrix(ht2_t))) %*% transmat)
          ht <- ht1+ht2
        }
        
        Bt[[c]] <- Bt[[c]] + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% solve(as.matrix(Matrix::bdiag(vc))) %*% as.matrix(diag(1/rep(as.vector(Wi),num_feature))) %*% as.matrix(Matrix::bdiag(dmuc)) - tau[i,c]*ht%o%qt
        
        #####################
        
        # Bg[[c]] <- Bg[[c]] + int_outer*tau[i,c] - tau[i,c]*qzeta%o%qzeta
        # Bt[[c]] <- Bt[[c]] + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% solve(as.matrix(Matrix::bdiag(vc))) %*% as.matrix(Matrix::bdiag(dmuc)) - tau[i,c]*qt%o%qt
        # Bgl[[c]] <- Bgl[[c]] - tau[i,c]*qzeta%o%qlambda
        
        Bag[[c]] <- Bag[[c]] - tau[i,c]*qalpha%o%qzeta
        Bal[[c]] <- Bal[[c]] - tau[i,c]*qalpha%o%qlambda
        Bat[[c]] <- Bat[[c]] - tau[i,c]*qalpha%o%qt
        Bta[[c]] <- Bta[[c]] - tau[i,c]*ht%o%qalpha
        
        if (length(mu[[c]][[j]]) == 1){
          Btg[[c]] <- Btg[[c]] - tau[i,c]*ht%o%qzeta
          Btl[[c]] <- Btl[[c]] - tau[i,c]*ht%o%qlambda
        }else{
          Btg[[c]] <- Btg[[c]] - tau[i,c]*ht%o%qzeta + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% (((as.vector(yi)-muc)*as.vector(dWi/Wi)) * do.call("rbind", rep(list(rbind(0,Hi[1:(nrow(yi)-1),])), num_feature ))) 
          Btl[[c]] <- Btl[[c]] - tau[i,c]*ht%o%qlambda + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% (((as.vector(yi)-muc)/as.vector(Wi)) * do.call("rbind", rep(list(longiterm), num_feature ))) 
        }
        
        Bgt[[c]] <- Bgt[[c]] - tau[i,c]*qzeta%o%qt
        Blt[[c]] <- Blt[[c]] - tau[i,c]*qlambda%o%qt
        
        bigBa = bigBa + tau[i,c]*qalpha
        bigBg[[c]] = tau[i,c]*qzeta
        bigBl[[c]] = tau[i,c]*qlambda
        bigBt[[c]] = tau[i,c]*qt
        bigBht[[c]] = tau[i,c]*ht
        
        
      }
      
      
      
    }
    
  }
  
  for (i in 1:n){
    # print(i)
    bigBa <- 0
    bigBg <- list()
    bigBl <- list()
    bigBt <- list()
    bigBht <- list()
    x = as.matrix(baseline[i,covx])
    xi = c(1,x)
    Hi = as.matrix(dat[dat$newid==i,c(covx,covy)])
    yi = as.matrix(dat[dat$newid==i,covy])
    time <- dat$time[dat$newid==i]
    time1 <- dat$time1[dat$newid==i]
    delta <- dat$delta[dat$newid==i]
    Ba <- 0
    
    for (c in 1:num_class){
      
      exphat <- exp(Hi%*%zeta[,c])
      lagexphat <- c(1,exphat[1:(length(exphat)-1)])
      difft <- timevec[hazvec[,c]>0,c]
      chazc <- chazvec[timevec[,c] %in% time1,c]
      hazc <- hazvec[timevec[,c] %in% time1,c]
      diffchazc <- c(chazc[1],diff(chazc))
      
      Bag[[c]] <- 0
      Bal[[c]] <- 0
      Bat[[c]] <- 0
      Bta[[c]] <- 0
      Bg[[c]] <- 0
      Bgl[[c]] <- 0
      Btg[[c]] <- 0
      Bl[[c]] <- 0
      Btl[[c]] <- 0
      Bt[[c]] <- 0
      Blg[[c]] <- 0
      Bgt[[c]] <- 0
      Blt[[c]] <- 0
      
      
      Dalpha <- -outer(p[i,-1],p[i,-1])
      diag(Dalpha) <- p[i,-1]*(1-p[i,-1])
      qalpha <- rep(0,(num_class-1)*(ncol(x)+1))
      for (j in 2:num_class){
        qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] = qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] - p[i,j]*xi
      }
      if (c > 1){
        qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] = qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] + xi
      }
      Ba = Ba + tau[i,c]*kronecker(Dalpha,outer(xi,xi)) - tau[i,c]*qalpha%o%qalpha
      
      # print(Ba)
      
      # survival part
      
      qlambda <- rep(0,length(difft))
      Bgl[[c]] = matrix(0,nrow=length(zeta[,c]),ncol=length(difft))
      Bl[[c]] = matrix(0,nrow=length(difft),ncol=length(difft)) #tau[i,c]*diag(1/(hazvec[hazvec[,c]>0,c]^2))
      int_scalar = 0
      int_vector = 0
      int_outer = 0
      
      qzeta <- as.vector(t(Hi) %*% (delta - exphat*diffchazc))
      
      qlambda <- as.vector(sapply(time1,function(x) x == difft) %*% ifelse(delta==1,1,0) - 
                             (sapply(time1,function(x) x >= difft)*sapply(time,function(x) x <= difft)*hazvec[timevec[,c] %in% difft,c]) %*% exphat)
      
      # qlambda <- as.vector(sapply(time1,function(x) x == difft) %*% ifelse(delta==1,1/hazc,0) - (sapply(time1,function(x) x >= difft)*sapply(time,function(x) x <= difft)) %*% exphat)
      
      
      Bl[[c]] = Bl[[c]] + tau[i,c]*diag(as.vector((sapply(time1,function(x) x >= difft)*sapply(time,function(x) x <= difft)) %*% exphat))
      for (j in 1:length(time1)){
        Bg[[c]] = Bg[[c]] + tau[i,c]*outer(Hi[j,],Hi[j,])*exphat[j]*diffchazc[j]
        Bgl[[c]] = Bgl[[c]] + tau[i,c]*exphat[j]*Hi[j,] %o% as.vector(sapply(time1[j],function(x) x >= difft)*sapply(time[j],function(x) x <= difft))
        Blg[[c]] = Blg[[c]] + tau[i,c]*exphat[j]*as.vector(sapply(time1[j],function(x) x >= difft)*sapply(time[j],function(x) x <= difft)*hazvec[timevec[,c] %in% difft,c]) %o% Hi[j,]
      }
      
      Bg[[c]] = Bg[[c]] - tau[i,c]*qzeta%o%qzeta
      Bl[[c]] = Bl[[c]] - tau[i,c]*qlambda%o%qlambda
      Bgl[[c]] = Bgl[[c]] - tau[i,c]*qzeta%o%qlambda
      Blg[[c]] = Blg[[c]] - tau[i,c]*qlambda%o%qzeta
      
      if (length(time) > 1){
        longiterm <- matrix(0,nrow=length(time),ncol=length(difft))
        longiterm = (t(sapply(time1,function(x) x >= difft)*sapply(time,function(x) x <= difft)) * as.vector(c(lagexphat[-1],0)))
      }
      
      # longitudinal part
      
      mu[[c]] <- list()
      v[[c]] <- list()
      dmu[[c]] <-list()
      dv[[c]] <- list()
      Wi = as.vector(W[dat$newid==i,c])
      dWi = as.vector(dW[dat$newid==i,c])
      
      for (j in 1:num_feature){
        if (Y_dist[j] == 'normal'){
          mulink <- function(x) x
          varlink <- function(x) rep(1,length(x))
          dmulink <- function(x) 1
          dvarlink <- function(x) 0
        }else if (Y_dist[j] == 'poi'){
          mulink <- function(x) exp(x)
          varlink <- function(x) x
          dmulink <- function(x) exp(x)
          dvarlink <- function(x) 1
        }else if (Y_dist[j] == 'bin'){
          mulink <- function(x) exp(x)/(1+exp(x))
          varlink <- function(x) x*(1-x)
          dmulink <- function(x) exp(x)/(1+exp(x))^2 #(exp(x)*(1+exp(x)) - exp(x)*exp(x))/(1+exp(x))^2
          dvarlink <- function(x) 1-2*x
        }
        
        mu[[c]][[j]] <- Mu[dat$newid==i,j,c]
        if (length(mu[[c]][[j]]) == 1){
          v[[c]][[j]] <- phi[j,c]*diag(sqrt(varlink(mu[[c]][[j]])),nrow=1) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(sqrt(varlink(mu[[c]][[j]])),nrow=1)
        }else{
          v[[c]][[j]] <- phi[j,c]*diag(as.vector(sqrt(varlink(mu[[c]][[j]])))) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j]]))))
        }
        dmu[[c]][[j]] <- cbind(dmulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$newid==i]))),
                               dmulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$newid==i])))*dat$time[dat$newid==i])
        
        if (length(mu[[c]][[j]]) == 1){
          dv_uv <- phi[j,c]*diag(dvarlink(mu[[c]][[j]])/(2*sqrt(varlink(mu[[c]][[j]]))),nrow=1) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(sqrt(varlink(mu[[c]][[j]])),nrow=1)
          dv_uv_dmu <- outer(dv_uv,dmu[[c]][[j]])
          dv[[c]][[j]] <- array(0,dim=c(length(mu[[c]][[j]]),length(mu[[c]][[j]]),2)) # here I only consider intercept and time as the covariates (2 dimensional)
          dv[[c]][[j]][1,1,] = (dv_uv_dmu[1,1,1,] + dv_uv_dmu[1,1,1,])
        }else{
          dv_uv <- phi[j,c]*diag(as.vector(dvarlink(mu[[c]][[j]])/(2*sqrt(varlink(mu[[c]][[j]]))))) %*% (gamma[j,c]^as.matrix(dist(dat$time[dat$newid==i]))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j]]))))
          dv_uv_dmu <- outer(dv_uv,dmu[[c]][[j]])
          dv[[c]][[j]] <- array(0,dim=c(length(mu[[c]][[j]]),length(mu[[c]][[j]]),2)) # here I only consider intercept and time as the covariates (2 dimensional)
          for (uu in 1:length(mu[[c]][[j]])){
            for (vv in 1:length(mu[[c]][[j]])){
              dv[[c]][[j]][uu,vv,] = (dv_uv_dmu[uu,vv,uu,] + dv_uv_dmu[vv,uu,vv,])
            }
          }
        }
      }
      
      muc <- unlist(mu[[c]])
      mu01 <- unlist(mu[[1]])
      dmuc <- dmu[[c]]
      vc <- v[[c]]
      dmu1 <- dmu[[1]]
      v1 <- v[[1]]
      
      dv1_int <- Matrix::bdiag(lapply(dv[[1]],function(x) (x[,,1])))
      dv1_t <- Matrix::bdiag(lapply(dv[[1]],function(x) (x[,,2])))
      dvc_int <- Matrix::bdiag(lapply(dv[[c]],function(x) (x[,,1])))
      dvc_t <- Matrix::bdiag(lapply(dv[[c]],function(x) (x[,,2])))
      vinvc <- as.matrix(solve(Matrix::bdiag(vc)))
      vinv1 <- as.matrix(solve(Matrix::bdiag(v1)))
      
      qt <- as.vector(t(as.matrix(Matrix::bdiag(dmuc))) %*% vinvc %*% ((as.vector(yi)-muc)/as.vector(Wi)) )
      
      if (c == 1){
        ht = rep(0,num_feature*2)
      }else{
        ht1 <- as.vector(t(as.matrix(Matrix::bdiag(dmuc)) - as.matrix(Matrix::bdiag(dmu1))) %*% (vinvc %*% ((as.vector(yi)-muc)) + vinv1 %*% ((as.vector(yi)-mu01))))/(2*log(ew[i,c]))
        ht2_int <- (mu01 - muc) * (- vinvc %*% dvc_int %*% vinvc %*% ((as.vector(yi)-muc))
                                   - vinv1 %*% dv1_int %*% vinv1 %*% ((as.vector(yi)-mu01)))/(2*log(ew[i,c]))
        ht2_t <- (mu01 - muc) * (- vinvc %*% dvc_t %*% vinvc %*% ((as.vector(yi)-muc))
                                 - vinv1 %*% dv1_t %*% vinv1 %*% ((as.vector(yi)-mu01)))/(2*log(ew[i,c]))
        transmat <- kronecker(diag(6),rep(1,length(mu[[c]][[j]])))
        # print(i)
        # print(ht2_int)
        # print(ht2_t)
        ht2 = as.vector(t(cbind(as.matrix(ht2_int),as.matrix(ht2_t))) %*% transmat)
        ht <- ht1+ht2
      }
      
      Bt[[c]] <- Bt[[c]] + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% solve(as.matrix(Matrix::bdiag(vc))) %*% as.matrix(diag(1/rep(as.vector(Wi),num_feature))) %*% as.matrix(Matrix::bdiag(dmuc)) - tau[i,c]*ht%o%qt
      
      #####################
      
      # Bg[[c]] <- Bg[[c]] + int_outer*tau[i,c] - tau[i,c]*qzeta%o%qzeta
      # Bt[[c]] <- Bt[[c]] + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% solve(as.matrix(Matrix::bdiag(vc))) %*% as.matrix(Matrix::bdiag(dmuc)) - tau[i,c]*qt%o%qt
      # Bgl[[c]] <- Bgl[[c]] - tau[i,c]*qzeta%o%qlambda
      
      Bag[[c]] <- Bag[[c]] - tau[i,c]*qalpha%o%qzeta
      Bal[[c]] <- Bal[[c]] - tau[i,c]*qalpha%o%qlambda
      Bat[[c]] <- Bat[[c]] - tau[i,c]*qalpha%o%qt
      Bta[[c]] <- Bta[[c]] - tau[i,c]*ht%o%qalpha
      
      if (length(mu[[c]][[j]]) == 1){
        Btg[[c]] <- Btg[[c]] - tau[i,c]*ht%o%qzeta
        Btl[[c]] <- Btl[[c]] - tau[i,c]*ht%o%qlambda
      }else{
        Btg[[c]] <- Btg[[c]] - tau[i,c]*ht%o%qzeta + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% (((as.vector(yi)-muc)*as.vector(dWi/Wi)) * do.call("rbind", rep(list(rbind(0,Hi[1:(nrow(yi)-1),])), num_feature ))) 
        Btl[[c]] <- Btl[[c]] - tau[i,c]*ht%o%qlambda + tau[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% (((as.vector(yi)-muc)/as.vector(Wi)) * do.call("rbind", rep(list(longiterm), num_feature ))) 
      }
      
      Bgt[[c]] <- Bgt[[c]] - tau[i,c]*qzeta%o%qt
      Blt[[c]] <- Blt[[c]] - tau[i,c]*qlambda%o%qt
      
      bigBa = bigBa + tau[i,c]*qalpha
      bigBg[[c]] = tau[i,c]*qzeta
      bigBl[[c]] = tau[i,c]*qlambda
      bigBt[[c]] = tau[i,c]*qt
      bigBht[[c]] = tau[i,c]*ht
      
    }
    
    BBa = BBa + Ba + outer(bigBa,bigBa)
    BBag = BBag + do.call("cbind",Bag) + outer(bigBa,unlist(bigBg))
    BBal = BBal + do.call("cbind",Bal) + outer(bigBa,unlist(bigBl))
    BBat = BBat + do.call("cbind",Bat) + outer(bigBa,unlist(bigBt))
    
    BBg = BBg + Matrix::bdiag(Bg) + outer(unlist(bigBg),unlist(bigBg)) # 
    BBl = BBl + Matrix::bdiag(Bl) + outer(unlist(bigBl),unlist(bigBl)) #
    BBgl = BBgl + Matrix::bdiag(Bgl) + outer(unlist(bigBg),unlist(bigBl))
    BBlg = BBlg + Matrix::bdiag(Blg) + outer(unlist(bigBl),unlist(bigBg))
    
    BBt = BBt + Matrix::bdiag(Bt) + outer(unlist(bigBht),unlist(bigBt))
    BBtg = BBtg + Matrix::bdiag(Btg) + outer(unlist(bigBht),unlist(bigBg))
    BBtl = BBtl + Matrix::bdiag(Btl) + outer(unlist(bigBht),unlist(bigBl))
    BBta = BBta + do.call("rbind",Bta) + outer(unlist(bigBht),bigBa)
    BBgt = BBgt + Matrix::bdiag(Bgt) + outer(unlist(bigBg),unlist(bigBt))
    BBlt = BBlt + Matrix::bdiag(Blt) + outer(unlist(bigBl),unlist(bigBt))
    
    S = S + outer(c(bigBa,unlist(bigBg),unlist(bigBl),unlist(bigBt)),c(bigBa,unlist(bigBg),unlist(bigBl),unlist(bigBt)))
  }
  
  BBa = as.matrix(BBa)
  BBag = as.matrix(BBag)
  BBal = as.matrix(BBal)
  BBat = as.matrix(BBat)
  
  BBg = as.matrix(BBg)
  BBl = as.matrix(BBl)
  BBgl = as.matrix(BBgl)
  BBlg = as.matrix(BBlg)
  
  BBt = as.matrix(BBt)
  BBtg = as.matrix(BBtg)
  BBtl = as.matrix(BBtl)
  BBgt = as.matrix(BBgt)
  BBlt = as.matrix(BBlt)
  BBta = as.matrix(BBta)
  
  # list(BBtg)
  
  I = rbind(cbind(BBa,BBag,BBal,BBat),
            cbind(t(BBag),BBg,BBgl,BBgt),
            cbind(t(BBal),BBlg,BBl,BBlt),
            cbind(BBta,BBtg,BBtl,BBt))
  
  # I = rbind(cbind(BBa,BBag,BBal,BBat),
  #           cbind(t(BBag),BBg,BBgl,t(BBtg)*0),
  #           cbind(t(BBal),t(BBgl),BBl,t(BBtl)*0),
  #           cbind(t(BBat),BBtg*0,BBtl*0,BBt))
  
  # I = rbind(cbind(BBa,BBag*0,BBal*0,BBat*0),
  #           cbind(t(BBag)*0,BBg,BBgl,t(BBtg)*0),
  #           cbind(t(BBal)*0,t(BBgl),BBl,t(BBtl)*0),
  #           cbind(t(BBat)*0,BBtg*0,BBtl*0,BBt))
  
  Sigma1 = solve(I+diag(ncol(I))*1e-4)
  Sigma2 = Sigma1%*%S%*%t(Sigma1)
  
  ASE <- sqrt(diag(Sigma2))
  ASEalpha_rb <- ASE[1:((num_class-1)*(ncol(x)+1))]
  ASEzeta_rb <- ASE[((num_class-1)*(ncol(x)+1)+1):((num_class-1)*(ncol(x)+1)+length(zeta))]
  ASEbeta_rb <- ASE[(length(ASE)-2*length(beta0)+1):length(ASE)]
  ASEbeta0_rb <- ASEbeta_rb[seq(1,num_feature*2*num_class,2)]
  ASEbeta1_rb <- ASEbeta_rb[seq(2,num_feature*2*num_class+1,2)]
  
  list(ASEalpha_rb,ASEzeta_rb,ASEbeta0_rb,ASEbeta1_rb)
  
}