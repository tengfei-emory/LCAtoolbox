// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

Rcpp::List loglik_EM(arma::mat alpha, arma::vec zeta, arma::mat x, arma::vec delta, arma::vec t, arma::vec tevent, double num_class, arma::vec d, arma::mat d1, arma::cube d2, arma::mat tau){
  int n = t.n_elem;
  
  arma::vec l_individual = zeros<vec>(n);
  arma::vec L_individual = zeros<vec>(n);
  double l = 0;
  double L = 0;
  arma::vec la = zeros<vec>(alpha.n_cols*alpha.n_rows);
  arma::vec lz = zeros<vec>(zeta.n_elem);
  //arma::vec lz1 = zeros<vec>(zeta.n_elem);
  
  arma::mat Ia = zeros<mat>(la.n_elem,la.n_elem);
  arma::mat Iz = zeros<mat>(lz.n_elem,lz.n_elem);
  //arma::mat Iz1 = zeros<mat>(lz.n_elem,lz.n_elem);
  arma::mat S = zeros<mat>(la.n_elem+lz.n_elem,la.n_elem+lz.n_elem);
  //arma::mat S1 = zeros<mat>(la.n_elem+lz.n_elem,la.n_elem+lz.n_elem);
  
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
    if (i == 0){
      idx = find(t < tevent(i+1));
    }else if (i < nt-1){
      arma::uvec lidx1 = find(t >= tevent(i)); 
      arma::uvec lidx2 = find(t < tevent(i+1)); 
      idx = intersect(lidx1,lidx2);
    }else{
      idx = find(t >= tevent(i));
    }
    int jidx = idx.n_elem;
    
    for (int j = 0; j < jidx; ++j){
      
      //arma::uvec id = idx(j);
      arma::vec xi = zeros<vec>(1+x.n_cols);
      xi(0) = 1;
      xi.subvec(1,x.n_cols) = x.row(idx(j)).t();
      arma::vec Sj = zeros<vec>(la.n_elem+lz.n_elem);
      //arma::vec Sj1 = zeros<vec>(la.n_elem+lz.n_elem);
      
      for (int c = 0; c < num_class; ++c){

        arma::vec zi = zeros<vec>(num_class-1+x.n_cols*num_class);
        zi.subvec(num_class-1,num_class-1+x.n_cols-1) = x.row(idx(j)).t();
        if (c > 0){
          zi(c-1) = 1;
          zi.subvec(num_class+c*x.n_cols-1,num_class+(c+1)*x.n_cols-1-1) = x.row(idx(j)).t();
        }
        
        arma::vec qalpha = zeros<vec>((num_class-1)*(x.n_cols+1));
        
        
        for (int k = 1; k < num_class; ++k){
          qalpha.subvec(k*(x.n_cols+1)-x.n_cols-1,k*(x.n_cols+1)-1) = - p(idx(j),k)*xi;
        }
        if (c > 0){
          qalpha.subvec(c*(x.n_cols+1)-x.n_cols-1,c*(x.n_cols+1)-1) = qalpha.subvec(c*(x.n_cols+1)-x.n_cols-1,c*(x.n_cols+1)-1) + xi;
        }
        
        arma::mat psub = p.row(idx(j));
        psub = psub.cols(1,p.n_cols-1);
        arma::mat Dalpha = -psub.t()*psub;
        Dalpha.diag() = psub % (1-psub);
        
        arma::mat pxi = xi*xi.t();
        Dalpha = kron(-Dalpha,pxi);
        
        double ds = sum(d.subvec(0,i));
        arma::mat d1s = sum(d1.rows(0,i),0);
        arma::mat d2s = sum(d2.slices(0,i),2);
        
        if (t(idx(j)) < tevent(i)){
          ds = 0;
          d1s = zeros<mat>(1,zeta.n_elem);
          d2s = zeros<mat>(zeta.n_elem,zeta.n_elem);
        }
        
        arma::vec qzeta1 = zi*delta(idx(j)) + delta(idx(j))*d1.row(i).t()/d(i);
        //arma::vec qzeta = qzeta1 - exp(dot(zi,zeta))*d1s.t() - zi*exp(dot(zi,zeta)) * ds;
        
        arma::mat Dzeta1 =  - delta(idx(j))*((d1.row(i).t()*d1.row(i))/pow(d(i),2) -  d2.slice(i)/d(i));
        //arma::mat Dzeta = Dzeta1 + exp(dot(zi,zeta))*(-(zi*zi.t())*ds - zi*d1s - d1s.t()*zi.t() - d2s);
        
        l_individual(idx(j)) = tau(idx(j),c)*( delta(idx(j)) * (log(d(i)) + dot(zi,zeta)) + log(p(idx(j),c)));
        
        if (t(idx(j)) < tevent(i)){
          l_individual(idx(j)) = tau(idx(j),c)*log(p(idx(j),c));
          L_individual(idx(j)) = L_individual(idx(j)) + p(idx(j),c);
        }else{
          L_individual(idx(j)) = L_individual(idx(j)) + p(idx(j),c)*pow(d(i)*exp(dot(zi,zeta)),delta(idx(j)))*exp(-ds*exp(dot(zi,zeta)));
        }
        
        //l0 = l0 + tau(idx(j),c)*( delta(idx(j)) * (log(d(i)) + dot(zi,zeta))  + log(p(idx(j),c)));;
        l = l + l_individual(idx(j));
        la = la + tau(idx(j),c)*qalpha;
        lz = lz + tau(idx(j),c)*qzeta1;
        //lz1 = lz1 + tau(idx(j),c)*qzeta;
        
        Ia = Ia + tau(idx(j),c)*Dalpha;
        Iz = Iz + tau(idx(j),c)*Dzeta1;
        //Iz1 = Iz1 + tau(idx(j),c)*Dzeta;
        
        arma::vec q = join_cols(qalpha,qzeta1);
        Sj = Sj + tau(idx(j),c)*q;
        //arma::vec q1 = join_cols(qalpha,qzeta);
        //Sj1 = Sj1 + tau(idx(j),c)*q1;
      }
      
      S = S + Sj*Sj.t();
      //S1 = S1 + Sj1*Sj1.t();
      
    }
    
  }
  
  L = sum(log(L_individual));
  
  arma::vec score = join_cols(la,lz);
  //arma::vec score1 = join_cols(la,lz1);
  
  arma::mat I = zeros<mat>(la.n_elem+lz.n_elem,la.n_elem+lz.n_elem);
  I.submat(0,0,la.n_elem-1,la.n_elem-1) = Ia;
  I.submat(la.n_elem,la.n_elem,I.n_rows-1,I.n_cols-1) = Iz;
  //Rcpp::Rcout << "flag1" << I;
  
  //arma::mat I1 = zeros<mat>(la.n_elem+lz.n_elem,la.n_elem+lz.n_elem);
  //I1.submat(0,0,la.n_elem-1,la.n_elem-1) = Ia;
  //I1.submat(la.n_elem,la.n_elem,I.n_rows-1,I.n_cols-1) = Iz1;
  
  Rcpp::List ret;
  ret["loglik"] = l;
  ret["obsloglik"] = L;
  //ret["loglik0"] = l0;
  ret["score"] = score;
  //ret["score1"] = score1;
  ret["I"] = I;
  //ret["I1"] = I1;
  ret["S"] = S;
  //ret["S1"] = S1;
  ret["scorez"] = lz;
  ret["Iz"] = Iz;
  ret["l_indiv"] = l_individual;
  ret["L_indiv"] = L_individual;
  return ret;
}
