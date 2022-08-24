// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

arma::vec newton(double maxiterNewton, double tolNewton, arma::mat alpha, arma::vec zeta, arma::mat x, arma::vec delta, arma::vec t, arma::vec tevent, double num_class, arma::mat tau, Rcpp::Function lambda_EM, Rcpp::Function loglik_EM){
  
  arma::vec zeta0 = zeta;
  double diff = 10;
  int count = 0;
  
  while (diff > tolNewton){
    count = count + 1;
    //Rcpp::Rcout << "count: " << count << endl;
    Rcpp::List lambda_res = lambda_EM(zeta,x,t,tevent,num_class,tau);
    arma::vec d = lambda_res["d"];
    arma::mat d1 = lambda_res["d1"];
    arma::cube d2 = lambda_res["d2"];
    Rcpp::List res = loglik_EM(alpha,zeta,x,delta,t,tevent,num_class,d,d1,d2,tau);
    arma::vec score = res["scorez"];
    arma::mat I = res["Iz"];
    arma::mat epsilon = arma::eye(I.n_rows,I.n_cols)*0.0001;
    arma::mat identmat = arma::eye(I.n_rows,I.n_cols);
    arma::vec step = solve(I+epsilon,identmat)*score;
    diff = max(abs(step));
    //Rcpp::Rcout << "count: " << count << step << diff << endl;
    if (diff > 10) break;
    zeta = zeta - step;
    //Rcpp::Rcout << "zeta: " << zeta << endl;
    if (count > maxiterNewton) break;
  }
  
  return zeta;
}