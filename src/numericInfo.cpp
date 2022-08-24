// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat numericInfo(arma::mat alpha, arma::vec zeta, arma::mat x, arma::vec delta, arma::vec t, arma::vec tevent, double num_class, arma::vec l0, arma::mat tau, double epsilon, Rcpp::Function lambda_EM, Rcpp::Function loglik_EM, Rcpp::Function postweight){

  Rcpp::List lambda_res = lambda_EM(zeta,x,t,tevent,num_class,tau);
  arma::vec d = lambda_res["d"];
  arma::mat d1 = lambda_res["d1"];
  arma::cube d2 = lambda_res["d2"];
  
  int n = t.n_elem;
  arma::vec alphavec = vectorise(alpha);
  arma::mat nInfo = zeros<mat>(alphavec.n_elem+zeta.n_elem,alphavec.n_elem+zeta.n_elem);
  arma::mat nlind = zeros<mat>(n,alphavec.n_elem+zeta.n_elem);
  //arma::vec nScore = zeros<vec>(alphavec.n_elem+zeta.n_elem);
  
  for (int i=0; i < alphavec.n_elem; ++i){
    
    arma::vec alphainc = vectorise(alpha);
    alphainc(i) = alphainc(i) + epsilon;
    arma::mat alphaincm = reshape(alphainc,x.n_cols+1,num_class-1);
    //Rcpp::List taures = postweight(alphaincm,zeta,x,delta,t,tevent,num_class,d);
    //arma::mat tau1 = taures["tau"];
    //Rcpp::List lambda_res_zinc = lambda_EM(zeta,x,t,tevent,num_class,tau1);
    //arma::vec dz = lambda_res_zinc["d"];
    //arma::mat d1z = lambda_res_zinc["d1"];
    //arma::cube d2z = lambda_res_zinc["d2"];
    Rcpp::List res_inc = loglik_EM(alphaincm,zeta,x,delta,t,tevent,num_class,d,d1,d2,tau);
    arma::vec resinclind = res_inc["l_indiv"];
    nlind.col(i) = resinclind;
    
  }
  

  
  for (int j=0; j < zeta.n_elem; ++j){
    
    arma::vec zetainc = zeta;
    
    zetainc(j) = zetainc(j) + epsilon;
    
    Rcpp::List lambda_res_zinc = lambda_EM(zetainc,x,t,tevent,num_class,tau);
    arma::vec dz = lambda_res_zinc["d"];
    arma::mat d1z = lambda_res_zinc["d1"];
    arma::cube d2z = lambda_res_zinc["d2"];
    //Rcpp::List taures = postweight(alpha,zetainc,x,delta,t,tevent,num_class,dz);
    //arma::mat tau1 = taures["tau"];
    //Rcpp::List lambda_res_zincnew = lambda_EM(zetainc,x,t,tevent,num_class,tau1);
    //arma::vec dznew = lambda_res_zincnew["d"];
    //arma::mat d1znew = lambda_res_zincnew["d1"];
    //arma::cube d2znew = lambda_res_zincnew["d2"];
    Rcpp::List res_inc = loglik_EM(alpha,zetainc,x,delta,t,tevent,num_class,dz,d1z,d2z,tau);
    arma::vec resinclind = res_inc["l_indiv"];
    //Rcpp::Rcout << "flag1" << resinclind << endl;
    nlind.col(j+alphavec.n_elem) = resinclind;
  
  }
  
  for (int k=0; k < nlind.n_cols; ++k){
    nlind.col(k) = (nlind.col(k) - l0)/epsilon;
  }
  
  for (int i=0; i < n; ++i){
    nInfo = nInfo + nlind.row(i).t() * nlind.row(i);
    //nScore = nScore + nlind.row(i).t();
  }
  
  //Rcpp::List ret;
  //ret["nInfo"] = nInfo;
  //ret["nlind"] = nlind;
  return nInfo;
}