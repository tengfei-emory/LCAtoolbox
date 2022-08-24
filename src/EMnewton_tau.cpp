// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
  
  #include "RcppArmadillo.h"
  // [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

Rcpp::List EMNewton(double maxiterEM, double tolEM, double maxiterNewton, double tolNewton, arma::mat alpha, arma::vec zeta, arma::mat x, arma::vec delta, arma::vec t, arma::vec tevent, double num_class, arma::vec d, arma::mat tau, Rcpp::Function lambda_EM, Rcpp::Function loglik_EM, Rcpp::Function postweight, Rcpp::Function printplot){
  
  Rcpp::List lambda_res;
  Rcpp::List res;
  // double alphadiff;
  // double zetadiff;
  double diff;
  arma::mat p;
  
  for (int i = 0; i < maxiterEM ; ++i){
    
    arma::mat alpha0 = vectorise(alpha);
    arma::vec zeta0 = zeta;
    arma::vec d0 = d;
    
    for (int j = 0; j < maxiterNewton; ++j){
      //arma::vec taudiff = vectorise(tau1) - vectorise(tau);
      //double diffEM = sum(vectorise(pow(taudiff,2)));
      //tau = tau1;
      
      lambda_res = lambda_EM(zeta,x,t,tevent,num_class,tau);
      arma::vec d00 = lambda_res["d"];
      d = d00;
      arma::mat d1 = lambda_res["d1"];
      arma::cube d2 = lambda_res["d2"];
      res = loglik_EM(alpha,zeta,x,delta,t,tevent,num_class,d,d1,d2,tau);
      arma::vec score = res["score"];
      arma::mat I = res["I"];
      //Rcpp::Rcout << "I: " << I << endl;
      //arma::mat epsilon = arma::eye(I.n_rows,I.n_cols)*0.0001;
      arma::mat identmat = arma::eye(I.n_rows,I.n_cols);
      //arma::vec step = inv(I+epsilon)*score;
      arma::vec step = solve(I,identmat)*score;
      
      alpha = alpha0 - step.subvec(0,alpha0.n_elem-1);
      alpha = reshape(alpha,x.n_cols+1,num_class-1);
      zeta = zeta - step.subvec(alpha0.n_elem,step.n_elem-1);
      
      Rcpp::List taures = postweight(alpha,zeta,x,delta,t,tevent,num_class,d);
      arma::mat tau1 = taures["tau"];
      arma::mat p1 = taures["p"];
      p = p1;
      tau = tau1;
      
      if (max(abs(step)) < tolNewton) break;
      
    }
    
    //Rcpp::Rcout << "alpha: " << alpha << endl;
    //Rcpp::Rcout << "zeta: " << zeta << endl;
    
    double l = res["obsloglik"];
    
    diff = sum(abs(vectorise(alpha) - alpha0)) + sum(abs(zeta - zeta0)) + sum(abs(d - d0));
    
    //double diffEM = max(abs(join_cols(alphadiff,zetadiff)));
    //Rcpp::Rcout << "iteration: " << i << " obsloglik: " << l << " diff: " << diff << endl;
    printplot(i,res,maxiterEM);
    if (diff < tolEM) break;
    
  }
  
  Rcpp::List ret;
  ret["alpha"] = alpha;
  ret["zeta"] = zeta;
  ret["lambda_res"] = lambda_res;
  ret["diff"] = diff;
  ret["tau"] = tau;
  ret["p"] = p;
  return ret;
}