// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

Rcpp::List varestcpp(arma::vec zeta, arma::vec delta, arma::vec chaz, arma::mat tau, arma::mat p, arma::mat x, arma::vec tevent){
  
  double num_class = tau.n_cols;
  double n = tau.n_rows;
  arma::vec dchaz = diff(chaz);
  arma::vec haz = zeros<vec>(dchaz.n_elem+1);
  haz(0) = chaz(0);
  haz.subvec(1,dchaz.n_elem) = dchaz;
  //haz(1,)
  //arma::vec haz = {chaz(0), dchaz};
  arma::vec unichaz = unique(chaz);
  unichaz = nonzeros(unichaz);
  arma::vec nonzerohaz = nonzeros(haz);
  
  arma::mat Ba = zeros<mat>((num_class-1)*(x.n_cols+1),(num_class-1)*(x.n_cols+1));
  arma::mat Bg = zeros<mat>(zeta.n_elem,zeta.n_elem);
  //arma::mat Bl = zeros<mat>(tevent.n_elem,tevent.n_elem);
  arma::mat Bl = diagmat((pow(nonzerohaz,-2)));
  arma::mat Bag = zeros<mat>((num_class-1)*(x.n_cols+1),zeta.n_elem);
  arma::mat Bal = zeros<mat>((num_class-1)*(x.n_cols+1),tevent.n_elem);
  arma::mat Bgl = zeros<mat>(zeta.n_elem,tevent.n_elem);
  arma::mat S = zeros<mat>((num_class-1)*(x.n_cols+1)+zeta.n_elem+tevent.n_elem,(num_class-1)*(x.n_cols+1)+zeta.n_elem+tevent.n_elem);
  arma::vec zi = zeros<vec>(num_class-1+x.n_cols*num_class);
  
  //Rcpp::Rcout << "flag1" << endl;
  
  for (int i=0; i < n; ++i){
    
    arma::vec bigBa = zeros<vec>((num_class-1)*(x.n_cols+1));
    arma::vec bigBg = zeros<vec>(zeta.n_elem);
    arma::vec bigBl = zeros<vec>(tevent.n_elem);
    
    arma::vec xi = zeros<vec>(1+x.n_cols);
    xi(0) = 1;
    xi.subvec(1,x.n_cols) = x.row(i).t();   
    
    //Rcpp::Rcout << "flag3" << i << endl;
    
    for (int c=0; c < num_class; ++c){
      
      zi = zeros<vec>(num_class-1+x.n_cols*num_class);
      zi.subvec(num_class-1,num_class-1+x.n_cols-1) = x.row(i).t();
      if (c > 0){
        zi(c-1) = 1;
        zi.subvec(num_class+c*x.n_cols-1,num_class+(c+1)*x.n_cols-1-1) = x.row(i).t();
      }
    
      //alpha
      arma::vec qalpha = zeros<vec>((num_class-1)*(x.n_cols+1));
      for (int k = 1; k < num_class; ++k){
        qalpha.subvec(k*(x.n_cols+1)-x.n_cols-1,k*(x.n_cols+1)-1) = - p(i,k)*xi;
      }
      if (c > 0){
        qalpha.subvec(c*(x.n_cols+1)-x.n_cols-1,c*(x.n_cols+1)-1) = qalpha.subvec(c*(x.n_cols+1)-x.n_cols-1,c*(x.n_cols+1)-1) + xi;
      }
      arma::mat psub = p.row(i);
      psub = psub.cols(1,p.n_cols-1);
      arma::mat Dalpha = -psub.t()*psub;
      Dalpha.diag() = psub % (1-psub);
      arma::mat pxi = xi*xi.t();
      Dalpha = kron(Dalpha,pxi);
      bigBa = bigBa + tau(i,c)*qalpha;
      
      //zeta
      arma::vec qzeta = zi*(delta(i) - exp(dot(zi,zeta))*chaz(i));
      bigBg = bigBg + tau(i,c)*qzeta;
      arma::mat Dzeta = (zi*zi.t())*exp(dot(zi,zeta))*chaz(i);
        
      //lambda
      
      //Rcpp::Rcout << "flag4" << i << endl;
      
      arma::vec vdelta1 = zeros<vec>(tevent.n_elem);
      arma::uvec idx1 = arma::find(unichaz == chaz(i));
      vdelta1(idx1) = delta(i)*ones<vec>(idx1.n_elem);

      arma::vec vdelta2 = zeros<vec>(tevent.n_elem);
      arma::uvec idx2 = arma::find(unichaz <= chaz(i));
      vdelta2(idx2) = ones<vec>(idx2.n_elem);
      
      arma::vec qlambda = -vdelta2*exp(dot(zi,zeta));
      if (haz(i) > 0){
        qlambda = qlambda + (1/haz(i))*vdelta1;
      }
      bigBl = bigBl + tau(i,c)*qlambda;
      
      //zeta lambda
      arma::mat Dzl = (zi*vdelta2.t())*exp(dot(zi,zeta));
      
      //contributions
      Ba = Ba + tau(i,c)*Dalpha - tau(i,c)*(qalpha*qalpha.t());
      Bag = Bag - tau(i,c)*(qalpha*qzeta.t());
      Bal = Bal - tau(i,c)*(qalpha*qlambda.t());
      Bg = Bg + tau(i,c)*Dzeta - tau(i,c)*(qzeta*qzeta.t());
      Bgl = Bgl + tau(i,c)*Dzl - tau(i,c)*(qzeta*qlambda.t());
      Bl = Bl - tau(i,c)*(qlambda*qlambda.t()); 
    }
    
    arma::vec q = join_cols(bigBa,bigBg,bigBl);
    S = S + q*q.t();
  
  }
  
  //Rcpp::Rcout << "flag2" << endl;
  
  arma::mat I = join_cols(join_rows(Ba,Bag,Bal),join_rows(Bag.t(),Bg,Bgl),join_rows(Bal.t(),Bgl.t(),Bl));
  I = I + S;
  
  arma::mat identmat = arma::eye(I.n_rows,I.n_cols);
  arma::mat Iinv = solve(I,identmat);
  arma::mat Irb = Iinv*S*Iinv;
  
  arma::vec ASErb = sqrt(Irb.diag());
  arma::vec ASErb_alpha = ASErb.subvec(0,(num_class-1)*(x.n_cols+1)-1);
  arma::vec ASErb_zeta = ASErb.subvec((num_class-1)*(x.n_cols+1),((num_class-1)*(x.n_cols+1)+zi.n_elem)-1);
  
  arma::vec ASEi = sqrt(Iinv.diag());
  arma::vec ASEi_alpha = ASEi.subvec(0,(num_class-1)*(x.n_cols+1)-1);
  arma::vec ASEi_zeta = ASEi.subvec((num_class-1)*(x.n_cols+1),((num_class-1)*(x.n_cols+1)+zi.n_elem)-1);
  
  Rcpp::List ret;
  // ret["ASErb"] = ASErb;
  // ret["ASEi"] = ASEi;
  //ret["I"] = I;
  //ret["S"] = S;
  ret["alpharb"] = ASErb_alpha;
  ret["zetarb"] = ASErb_zeta;
  ret["alphai"] = ASEi_alpha;
  ret["zetai"] = ASEi_zeta;
  ret["I"] = Irb;
  return ret;
}