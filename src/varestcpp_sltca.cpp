// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include "sltca_common.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

Rcpp::List varestcpp_sltca(arma::vec t0, arma::mat x, arma::mat y,
                           arma::mat tau, arma::mat p,
                           Rcpp::StringVector Y_dist, arma::vec id, arma::vec obs,
                           arma::cube mu, double lbeta, arma::mat covgee, arma::mat phi,
                           arma::mat gamma, std::string cor){

  double n = tau.n_rows;
  double num_class = tau.n_cols;
  double num_feature = y.n_cols;

  arma::mat Ba = zeros<mat>((num_class-1)*(x.n_cols+1),(num_class-1)*(x.n_cols+1));
  arma::mat Bb = zeros<mat>(num_class*num_feature*lbeta,num_class*num_feature*lbeta);
  arma::mat Bab = zeros<mat>((num_class-1)*(x.n_cols+1),num_class*num_feature*lbeta);
  arma::mat S = zeros<mat>((num_class-1)*(x.n_cols+1)+num_class*num_feature*lbeta,
                           (num_class-1)*(x.n_cols+1)+num_class*num_feature*lbeta);

  for (int i = 0; i < n; ++i){

    arma::vec bigBa = zeros<vec>((num_class-1)*(x.n_cols+1));
    arma::vec bigBb = zeros<vec>(num_feature*num_class*lbeta);

    arma::vec xi = zeros<vec>(1+x.n_cols);
    xi(0) = 1;
    xi.subvec(1,x.n_cols) = x.row(i).t();
    arma::vec geex = covgee.row(i).t();

    arma::uvec ididx = find(id == i+1);
    arma::mat yj = y.rows(ididx);
    arma::vec t0j = t0(ididx);
    arma::mat distime = sltca_distime(obs(ididx));

    for (int c = 0; c < num_class; ++c){

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

      Ba = Ba + tau(i,c)*Dalpha - tau(i,c)*(qalpha*qalpha.t());
      bigBa = bigBa + tau(i,c)*qalpha;

      // GEE
      arma::mat mujl = mu.slice(c).rows(ididx);
      arma::vec qbeta = zeros<vec>(num_feature*lbeta);
      arma::mat Dbeta = zeros<mat>(num_feature*lbeta,num_feature*lbeta);

      for (int l = 0; l < num_feature; ++l){

        arma::vec qseg;
        arma::mat Dblock;
        sltca_qscore(yj.col(l), mujl.col(l), t0j, geex,
                     Rcpp::as<std::string>(Y_dist(l)), phi(l,c), gamma(l,c),
                     cor, distime, lbeta, qseg, Dblock);

        qbeta.subvec(lbeta*l,lbeta*(l+1)-1) = qseg;
        Dbeta.submat(lbeta*l,lbeta*l,lbeta*(l+1)-1,lbeta*(l+1)-1) = Dblock;
      }

      Bb.submat(c*num_feature*lbeta,c*num_feature*lbeta,(c+1)*num_feature*lbeta-1,(c+1)*num_feature*lbeta-1) = Bb.submat(c*num_feature*lbeta,c*num_feature*lbeta,(c+1)*num_feature*lbeta-1,(c+1)*num_feature*lbeta-1) + tau(i,c)*Dbeta - tau(i,c)*(qbeta*qbeta.t());
      bigBb.subvec(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) = bigBb.subvec(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) + tau(i,c)*qbeta;
      Bab.cols(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) = Bab.cols(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) - tau(i,c)*(qalpha*qbeta.t());
    }

    arma::vec q = join_cols(bigBa,bigBb);
    S = S + q*q.t();
  }

  arma::mat I = join_cols(join_rows(Ba,Bab),
                          join_rows(Bab.t(),Bb));
  I = I + S;

  arma::mat identmat = arma::eye(I.n_rows,I.n_cols);
  arma::mat Iinv = solve(I,identmat);
  arma::mat Irb = Iinv*S*Iinv.t();

  arma::vec ASErb = sqrt(Irb.diag());
  arma::vec ASErb_alpha = ASErb.subvec(0,(num_class-1)*(x.n_cols+1)-1);
  arma::vec ASErb_beta = ASErb.subvec(ASErb.n_rows-num_feature*num_class*lbeta,ASErb.n_rows-1);

  arma::vec ASEi = sqrt(Iinv.diag());
  arma::vec ASEi_alpha = ASEi.subvec(0,(num_class-1)*(x.n_cols+1)-1);
  arma::vec ASEi_beta = ASEi.subvec(ASEi.n_rows-num_feature*num_class*lbeta,ASEi.n_rows-1);

  Rcpp::List ret;
  ret["alpharb"] = ASErb_alpha;
  ret["betarb"] = ASErb_beta;
  ret["alphai"] = ASEi_alpha;
  ret["betai"] = ASEi_beta;
  return ret;
}
