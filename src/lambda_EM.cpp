// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
// [[Rcpp::export]]

Rcpp::List lambda_EM(arma::vec zeta, arma::mat x, arma::vec t, arma::vec tevent, double num_class, arma::mat tau){
  int n = t.n_elem;
  arma::vec d = zeros<vec>(tevent.n_elem);
  arma::mat d1 = zeros<mat>(tevent.n_elem,zeta.n_elem);
  arma::cube d2 = zeros<cube>(zeta.n_elem,zeta.n_elem,tevent.n_elem);
  
  arma::vec dall = zeros<vec>(n);
  arma::mat d1all = zeros<mat>(n,zeta.n_elem);
  arma::cube d2all = zeros<cube>(zeta.n_elem,zeta.n_elem,n);
  
  
  for (int j = 0; j < n; ++j){
    
    for (int c = 0; c < num_class; ++c){
      
      arma::vec zi = zeros<vec>(num_class-1+x.n_cols*num_class);
      
      zi.subvec(num_class-1,num_class-1+x.n_cols-1) = x.row(j).t();
      
      if (c > 0){
        zi(c-1) = 1;
        zi.subvec(num_class+c*x.n_cols-1,num_class+(c+1)*x.n_cols-1-1) = x.row(j).t();
      }
      

      dall(j) = dall(j) + tau(j,c) * exp(dot(zi,zeta));
      d1all.row(j) = d1all.row(j) + zi.t() * (tau.at(j,c) * exp(dot(zi,zeta)));
      d2all.slice(j) = d2all.slice(j) + (zi*zi.t()) * (tau.at(j,c) * exp(dot(zi,zeta)));
      
    }
    
  }
  
  int nt = tevent.n_elem;
  for (int i = 0; i < nt; ++i){
    
    arma::uvec idx = arma::find(t >= tevent(i));
    
    double ds = sum(dall(idx));
    arma::mat d1s = sum(d1all.rows(idx),0);
    arma::mat d2s = sum(d2all.slices(idx),2);
    
    d2.slice(i) = -d2s/pow(ds,2) + 2*(d1s.t()*d1s)/pow(ds,3);
    d1.row(i) = -d1s/pow(ds,2);
    d(i) = 1/ds;
    
  }
  
  Rcpp::List ret;
  ret["d"] = d;
  ret["d1"] = d1;
  ret["d2"] = d2;
  return ret;
}