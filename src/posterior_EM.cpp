// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
// [[Rcpp::export]]

Rcpp::List postweight(arma::mat alpha, arma::vec zeta, arma::mat x, arma::vec delta, arma::vec t, arma::vec tevent, double num_class, arma::vec d){
  int n = t.n_elem;
  
  arma::mat survlik = ones<mat>(n,num_class);
  arma::mat tau = zeros<mat>(n,num_class);
  
  arma::mat X = zeros<mat>(n,1+x.n_cols);
  X.col(0) = ones<vec>(n);
  X.cols(1,x.n_cols) = x;
  arma::mat p = zeros<mat>(n,num_class);
  p.col(0) = ones<vec>(n);
  
  for (int c = 1; c < num_class; ++c){
    p.col(c) = exp(X*alpha.col(c-1));
  }
  for (int i = 0; i < n; ++i){
    double psum = sum(p.row(i));
    p.row(i) = p.row(i)/psum;
  }
  
  for (int i = 0; i < p.n_rows; ++i){
    for (int j = 0; j < p.n_cols; ++j){
      if(p(i,j) < 1e-8){
        p(i,j) = 1e-8;
      }
    }
  }
  
  int nt = tevent.n_elem;
  for (int i = 0; i < nt; ++i){
    arma::uvec idx;
    if (i < nt-1){
      arma::uvec lidx1 = find(t >= tevent(i)); 
      arma::uvec lidx2 = find(t < tevent(i+1)); 
      idx = intersect(lidx1,lidx2);
    }else{
      idx = find(t >= tevent(i));
    }
    int jidx = idx.n_elem;
    
    for (int j = 0; j < jidx; ++j){
      
      arma::vec xi = zeros<vec>(1+x.n_cols);
      xi(0) = 1;
      xi.subvec(1,x.n_cols) = x.row(idx(j)).t();
      
      for (int c = 0; c < num_class; ++c){
        
        arma::vec zi = zeros<vec>(num_class-1+x.n_cols*num_class);
        zi.subvec(num_class-1,num_class-1+x.n_cols-1) = x.row(idx(j)).t();
        if (c > 0){
          zi(c-1) = 1;
          zi.subvec(num_class+c*x.n_cols-1,num_class+(c+1)*x.n_cols-1-1) = x.row(idx(j)).t();
        }
        double ds = sum(d.subvec(0,i));
        survlik(idx(j),c) = pow(exp(dot(zi,zeta)),delta(idx(j)))*exp(-exp(dot(zi,zeta))*ds);
        //Rcpp::Rcout << "flag1 " << pow(exp(dot(zi,zeta)),delta(idx(j))) << " delta " << delta(idx(j)) << endl;
      }

    }
    
  }
  
  arma::mat pew = p % survlik;
  for (int i = 0; i < n; ++i){
    double pewsum = sum(pew.row(i));
    tau.row(i) = pew.row(i)/pewsum;
  }
  
  //Rcpp::Rcout << "flag1" << tau.n_rows << tau.n_cols;
  
  for (int i = 0; i < tau.n_rows; ++i){
    for (int j = 0; j < tau.n_cols; ++j){
      if(tau(i,j) < 1e-8){
        tau(i,j) = 1e-8;
      }
    }
  }
  
  Rcpp::List ret;
  ret["tau"] = tau;
  ret["p"] = p;
  return ret;
}