ic <- function(geefit,Y_dist,yy){
  
  #weights(geefit) gives posterior membership probabilities tau
  #geefit$geese$gamma is the estimated scale parameter
  
  if(Y_dist == 'normal'){
    qlik = -0.5*(yy-fitted(geefit))*(yy-fitted(geefit))/geefit$geese$gamma - 0.5*log(geefit$geese$gamma)
    # v = geefit$geese$gamma 
  }else if (Y_dist == 'bin'){
    qlik = (yy*log(fitted(geefit)/(1-fitted(geefit))) + log(1 - fitted(geefit)))/geefit$geese$gamma - 0.5*log(geefit$geese$gamma*fitted(geefit)*(1-fitted(geefit)))
    # v = geefit$geese$gamma*fitted(geefit)*(1-fitted(geefit))
  }else if (Y_dist == 'poi'){
    qlik = (yy*log(fitted(geefit)) - fitted(geefit))/geefit$geese$gamma - 0.5*log(geefit$geese$gamma*fitted(geefit))
    # v = geefit$geese$gamma*fitted(geefit)
  }
  
  # eqlik = qlik - 0.5*sum(weights(geefit)*log(v))
  
  qlik
  
}