brier <- function(s=0,t,traindat,testdat,coxfit,lcfit){
  
  testdat <- testdat[order(testdat$tildet),]
  testdat = testdat[testdat$tildet >= s,]
  
  ## IPCW
  
  invGfit <- survfit(Surv(tildet,1-delta)~1,data=traindat)
  
  if (max(invGfit$time) < t){
    Gt <- invGfit$surv[which(invGfit$time == max(invGfit$time))]
  }else{
    Gt <- invGfit$surv[min(which(invGfit$time>=t))]
  }
  
  if (s==0){
    Gs=1
  }else{
    Gs <- invGfit$surv[min(which(invGfit$time>=s))]
  }
  G <- rep(NA,nrow(testdat))
  for (i in 1:length(G)){
    G[i] <- invGfit$surv[max(which(invGfit$time<=testdat$tildet[i]))]
  }
  G[is.na(G)] = 1
  W <- ifelse(testdat$tildet > t,1/(Gt/Gs),testdat$delta/(G/Gs))
  
  ## Cox
  
  y <- traindat$tildet[traindat$delta==1]
  h0 <- rep(NA, length(y))
  for(l in 1:length(y))
  {
    h0[l] <- 1 / sum(exp(as.matrix(traindat[traindat$tildet >= y[l],covx]) %*% coxfit$coefficients))
  }
  chaz <- cumsum(h0)
  coef <- as.numeric(exp(as.matrix(testdat[,covx]) %*% coxfit$coefficients))
  
  if (min(y) > t){
    Shatt = rep(1,nrow(testdat))
  }else{
    Shatt <- exp(-chaz[max(which(y <= t))]*coef)
  }
  
  if (min(y) > s){
    Shats = 1
  }else{
    Shats <- exp(-chaz[max(which(y <= s))]*coef)
  }
  Shat <- Shatt/Shats
  
  ST <- rep(NA,nrow(testdat))
  for (i in 1:length(G)){
    if (min(y) > testdat$tildet[i]){
      ST[i] = 1
    }else{
      ST[i] <- exp(-chaz[max(which(y <= testdat$tildet[i]))]*coef[i])
    }
  }
  ST = ST/Shats
  
  BS2_cox=rep(0,nrow(testdat))
  for (i in 1:nrow(testdat)){
    if (testdat$tildet[i] > t){
      BS2_cox[i] = (1-Shat[i])^2
    }else if(testdat$tildet[i] <= t & testdat$delta[i] == 1){
      BS2_cox[i] = (0-Shat[i])^2
    }else{
      BS2_cox[i] = (1-Shat[i])^2*(Shat[i]/ST[i]) + (0-Shat[i])^2*(1-Shat[i]/ST[i])
    }
  }
  
  BS1_cox = mean(W*(I(testdat$tildet > t)-Shat)^2)
  BS2_cox = mean(BS2_cox)
  
  ## Latent (2 classes)
  
  pred.lsca <- predict.LSCA(lcfit,testdat,timepoint=NULL)
  pred.lsca.t <- predict.LSCA(lcfit,testdat,timepoint=t)
  if (min(y) > s){
    predsurvp <- pred.lsca.t$s
    STp <- pred.lsca$s
  }else{
    pred.lsca.s <- predict.LSCA(lcfit,testdat,timepoint=s)
    predsurvp <- pred.lsca.t$s/pred.lsca.s$s
    STp <- pred.lsca$s/pred.lsca.s$s
  }
  
  BS2_latent=rep(0,nrow(testdat))
  for (i in 1:nrow(testdat)){
    if (testdat$tildet[i] > t){
      BS2_latent[i] = (1-predsurvp[i])^2
    }else if(testdat$tildet[i] <= t & testdat$delta[i] == 1){
      BS2_latent[i] = (0-predsurvp[i])^2
    }else{
      BS2_latent[i] = (1-predsurvp[i])^2*(predsurvp[i]/STp[i]) + (0-predsurvp[i])^2*(1-predsurvp[i]/STp[i])
    }
  }
  
  BS1_latent = mean(W*(I(testdat$tildet > t)-predsurvp)^2)
  BS2_latent = mean(BS2_latent) 
  
  output <- list()
  output$BS1_cox = BS1_cox  
  output$BS1_latent = BS1_latent
  
  output$BS2_cox = BS2_cox
  output$BS2_latent = BS2_latent
  
  output$IBS1_cox = BS1_cox*(Gs/Gt)
  output$IBS1_latent = BS1_latent*(Gs/Gt)
  
  output$IBS2_cox = BS2_cox*(Gs/Gt)
  output$IBS2_latent = BS2_latent*(Gs/Gt)
  
  output$IBS_weight = Gs/Gt
  return(output)
}
