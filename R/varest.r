VarEst <- function(zeta,delta,chaz,tau,p,dat,x){
  
  require(progress)
  #require(Matrix)
  num_class = ncol(tau)
  #nmax <- max(dat$num_obs)
  n <- length(unique(dat$id))
  haz <- c(chaz[1],diff(chaz))
  
  num_jumps = sum(haz>0)
  
  Ba <- 0
  Bag <- 0
  Bal <- 0
  Bg <- 0
  Bgl <- 0
  Bl = diag(1/(haz[haz>0]^2))
  S <- 0
  print('Variance estimation started:')
  pb <- progress_bar$new(total=n)
  
  for(i in 1:n){
    
    pb$tick()
    Sys.sleep(1/n)
    
    bigBa = 0
    bigBg = 0
    bigBl = 0
    xi <- c(1,x[i,])
    
    #qalpha =qalpha - kronecker(p[i,],xi)
    for (c in 1:num_class){
      
      zi <- c(rep(0,num_class-1), x[i,], rep(0,length(x[i,])*(num_class-1)))
      
      if (c > 1){
        zi[c-1] = 1
        zi[(num_class+(c-1)*length(x[i,])):(num_class+c*length(x[i,])-1)] <- x[i,]
      }
      
      #bigBi[[c]] <- (tau[i,c] * q)#%o%(tau[i,c] * q)
      #bigA[[c]] <- bigA[[c]] + p[i,c]*t(as.matrix(bdiag(dmuc))) %*% as.matrix(solve(bdiag(vc))) %*% as.matrix(bdiag(dmuc)) - tau[i,c]*q%o%q
      
      Dalpha <- -outer(p[i,-1],p[i,-1])
      diag(Dalpha) <- p[i,-1]*(1-p[i,-1])
      qalpha <- rep(0,(num_class-1)*(ncol(x)+1))
      for (j in 2:num_class){
        qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] = qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] - p[i,j]*c(1,x[i,])
      }
      if (c > 1){
        qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] = qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] + c(1,x[i,])
      }
      
      qzeta <- zi*(delta[i] - as.numeric(exp(zi %*% zeta))*chaz[i])
      
      vdelta <- rep(0,num_jumps)
      idx <- which(unique(chaz[chaz>0]) == chaz[i])
      vdelta[idx] = delta[i]
      if (haz[i] == 0){
        qlambda <- -as.numeric(exp(zi %*% zeta))*I(unique(chaz[chaz>0]) <= chaz[i])
      }else{
        qlambda <- (1/haz[i])*vdelta-as.numeric(exp(zi %*% zeta))*I(unique(chaz[chaz>0]) <= chaz[i])
      }
      
      bigBa = bigBa + tau[i,c]*qalpha
      bigBg = bigBg + tau[i,c]*qzeta
      bigBl = bigBl + tau[i,c]*qlambda
      
      Ba = Ba + tau[i,c]*kronecker(Dalpha,outer(xi,xi)) - tau[i,c]*qalpha%o%qalpha
      #Ba = Ba + p[i,c]*kronecker(Dalpha,outer(xi,xi)) - tau[i,c]*qalpha%o%qalpha
      
      Bag = Bag - tau[i,c]*qalpha%o%qzeta
      Bal = Bal - tau[i,c]*qalpha%o%qlambda
      
      Bg = Bg + outer(zi,zi)*tau[i,c]*as.numeric(exp(zi %*% zeta))*chaz[i] - tau[i,c]*qzeta%o%qzeta
      #Bg = Bg + outer(zi,zi)*p[i,c]*as.numeric(exp(zi %*% zeta))*chaz[i] - tau[i,c]*qzeta%o%qzeta
      
      Bgl = Bgl + tau[i,c]*as.numeric(exp(zi %*% zeta))*zi%o%I(unique(chaz[chaz>0]) <= chaz[i]) - tau[i,c]*qzeta%o%qlambda
      
      Bl = Bl - tau[i,c]*qlambda%o%qlambda
      
    }
    
    Ba = Ba + outer(bigBa,bigBa)
    Bag = Bag + outer(bigBa,bigBg)
    Bal = Bal + outer(bigBa,bigBl)
    
    Bg = Bg + outer(bigBg,bigBg)
    Bgl = Bgl + outer(bigBg,bigBl)
    
    Bl = Bl + outer(bigBl,bigBl)
    
    S = S + outer(c(bigBa,bigBg,bigBl),c(bigBa,bigBg,bigBl))
  }
  
  I = rbind(cbind(Ba,Bag,Bal),cbind(t(Bag),Bg,Bgl),cbind(t(Bal),t(Bgl),Bl))
  
  Sigma1 = solve(I+diag(ncol(I))*1e-4)
  Sigma2 = Sigma1%*%S%*%Sigma1
  #Sigmaa = solve(Ba+1e-4)%*%S[1:((num_class-1)*(ncol(x)+1)),1:((num_class-1)*(ncol(x)+1))]%*%solve(Ba)
  #Sigmag = solve(Bg+1e-4)%*%S[((num_class-1)*(ncol(x)+1)+1):((num_class-1)*(ncol(x)+1)+length(zi)),((num_class-1)*(ncol(x)+1)+1):((num_class-1)*(ncol(x)+1)+length(zi))]%*%solve(Bg)
  
  ASE <- sqrt(diag(Sigma2))
  ASEalpha_rb <- ASE[1:((num_class-1)*(ncol(x)+1))]
  ASEzeta_rb <- ASE[((num_class-1)*(ncol(x)+1)+1):((num_class-1)*(ncol(x)+1)+length(zi))]
  #ASElambda <- ASE[]
  
  ASE <- sqrt(diag(Sigma1[1:((num_class-1)*(ncol(x)+1)+length(zi)),1:((num_class-1)*(ncol(x)+1)+length(zi))]))
  ASEalpha_i <- ASE[1:((num_class-1)*(ncol(x)+1))]
  ASEzeta_i <- ASE[((num_class-1)*(ncol(x)+1)+1):((num_class-1)*(ncol(x)+1)+length(zi))]
  
  ASE <- sqrt(diag(solve(S+diag(ncol(S))*1e-4)))
  ASEalpha_s <- ASE[1:((num_class-1)*(ncol(x)+1))]
  ASEzeta_s <- ASE[((num_class-1)*(ncol(x)+1)+1):((num_class-1)*(ncol(x)+1)+length(zi))]
  
  #ASEalpha_nv <- sqrt(diag(Sigmaa))
  #ASEzeta_nv <- sqrt(diag(Sigmag))
  
  #ASEalpha_nvi <- sqrt(diag(solve(Ba)))
  #ASEzeta_nvi <- sqrt(diag(solve(Bg)))
  
  list(ASEalpha_rb,ASEzeta_rb,ASEalpha_i,ASEzeta_i,ASEalpha_s,ASEzeta_s)
  
}