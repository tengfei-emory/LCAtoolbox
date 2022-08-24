// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

Rcpp::List varestcpp_sltca(arma::vec t0, arma::mat x, arma::mat y,
                           arma::mat tau, arma::mat p, 
                           Rcpp::StringVector Y_dist, arma::vec id,
                           arma::cube mu, double lbeta, arma::mat covgee, arma::mat phi, 
                           arma::mat gamma){
  
  // arma::vec delta, arma::vec chaz, arma::mat tau, arma::mat p, 
  // arma::mat x, arma::mat bsp, arma::vec tevent, arma::vec d){
  
  // zeta,phi,gamma,Mu,t0,t1,tevent,delta,x,xy,y,bsp,Y_dist,ew,d,tau,p,id,last
  
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
    arma::mat distime = zeros<mat>(yj.n_rows,yj.n_rows);
    int numobs = yj.n_rows;
    arma::vec t0j = t0(ididx);
    
    if (numobs > 1){
      // arma::vec distvec = arma::diff(t0j);
      arma::vec distvec = arma::ones(numobs-1);
      // Rcpp::Rcout << distvec << endl;
      for (int jj = 0; jj < numobs-1 ; ++jj){
        // distime.submat(j,j+1,j,numobs-1) = arma::diff(timei,j+1);
        distime(jj,span(jj+1,numobs-1)) = cumsum(distvec.subvec(jj,numobs-2)).t();
        distime(span(jj+1,numobs-1),jj) = cumsum(distvec.subvec(jj,numobs-2));
      }
    }
    
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
      
      arma::mat muc = mu.slice(c);
      arma::mat mujl = muc.rows(ididx);
      arma::vec qbeta = zeros<vec>(num_feature*lbeta);
      arma::mat Dbeta = zeros<mat>(num_feature*lbeta,num_feature*lbeta);
      
      for (int l = 0; l < num_feature; ++l){
        
        arma::vec yjl = yj.col(l);
        arma::vec mujlc = mujl.col(l);
        
        arma::mat vjlc = zeros<mat>(numobs,numobs);
        arma::mat vlinkjlc = zeros<mat>(numobs,numobs);
        arma::mat corrjlc = zeros<mat>(numobs,numobs);
        
        arma::mat dmujlc = zeros<mat>(numobs,lbeta);
        
        // Rcpp::Rcout << mujlc << endl;
        
        if (Y_dist(l) == "normal"){
          for (int ll=0; ll < numobs; ++ll){
            vlinkjlc(ll,ll) = 1;
            for (int kk=0; kk < numobs; ++kk){
              corrjlc(ll,kk) = pow(gamma(l,c),distime(ll,kk));
            }
            
            dmujlc(ll,0) = 1;
            dmujlc(ll,1) = t0j(ll);
            for (int geeidx = 2; geeidx < lbeta; ++geeidx){
              dmujlc(ll,geeidx) = geex(geeidx-2);
            }
          }
          vjlc = phi(l,c)*vlinkjlc*corrjlc*vlinkjlc;
          // Rcpp::Rcout << dmujlc << endl;
          // Rcpp::Rcout << vjlc << endl;
        }else if (Y_dist(l) == "poi"){
          for (int ll=0; ll < numobs; ++ll){
            vlinkjlc(ll,ll) = sqrt(mujlc(ll));
            for (int kk=0; kk < numobs; ++kk){
              corrjlc(ll,kk) = pow(gamma(l,c),distime(ll,kk));
            }
            
            dmujlc(ll,0) = mujlc(ll);
            dmujlc(ll,1) = mujlc(ll)*t0j(ll);
            for (int geeidx = 2; geeidx < lbeta; ++geeidx){
              dmujlc(ll,geeidx) = mujlc(ll)*geex(geeidx-2);
            }
          }
          vjlc = phi(l,c)*vlinkjlc*corrjlc*vlinkjlc;
          
        }else if (Y_dist(l) == "bin"){
          for (int ll=0; ll < numobs; ++ll){
            vlinkjlc(ll,ll) = sqrt(mujlc(ll)*(1-mujlc(ll)));
            for (int kk=0; kk < numobs; ++kk){
              corrjlc(ll,kk) = pow(gamma(l,c),distime(ll,kk));
            }
            
            dmujlc(ll,0) = mujlc(ll)*(1-mujlc(ll));
            dmujlc(ll,1) = mujlc(ll)*(1-mujlc(ll))*t0j(ll);
            for (int geeidx = 2; geeidx < lbeta; ++geeidx){
              dmujlc(ll,geeidx) = mujlc(ll)*(1-mujlc(ll))*geex(geeidx-2);
            }
          }
          vjlc = phi(l,c)*vlinkjlc*corrjlc*vlinkjlc;
          
        }
        arma::mat identmat = arma::eye(vjlc.n_rows,vjlc.n_cols);
        qbeta.subvec(lbeta*l,lbeta*(l+1)-1) = dmujlc.t()*solve(vjlc,identmat)*(yjl-mujlc);
        Dbeta.submat(lbeta*l,lbeta*l,lbeta*(l+1)-1,lbeta*(l+1)-1) = dmujlc.t()*solve(vjlc,identmat)*dmujlc;
        
      }
      
      Bb.submat(c*num_feature*lbeta,c*num_feature*lbeta,(c+1)*num_feature*lbeta-1,(c+1)*num_feature*lbeta-1) = Bb.submat(c*num_feature*lbeta,c*num_feature*lbeta,(c+1)*num_feature*lbeta-1,(c+1)*num_feature*lbeta-1) + tau(i,c)*Dbeta - tau(i,c)*(qbeta*qbeta.t());
      // Rcpp::Rcout << "flag4 " << endl;
      
      bigBb.subvec(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) = bigBb.subvec(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) + tau(i,c)*qbeta;
      // Rcpp::Rcout << "flag5 " << endl;
      
      Bab.cols(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) = Bab.cols(c*num_feature*lbeta,(c+1)*num_feature*lbeta-1) - tau(i,c)*(qalpha*qbeta.t());
      // Rcpp::Rcout << "flag8 " << endl;
      
      // Rcpp::Rcout << "id" << i << "tau" << tau(i,c) << endl;
      // 
      // Rcpp::Rcout << "id" << i << "w" << Wjc << endl;
      // 
      // Rcpp::Rcout << "a" << tau(i,c)*qalpha << endl;
      // //
      // Rcpp::Rcout << "z" << tau(i,c)*qzeta << endl;
      // //
      // Rcpp::Rcout << "b" << tau(i,c)*qbeta << endl;
    }
    
    
    arma::vec q = join_cols(bigBa,bigBb);
    // qall = qall + q;
    S = S + q*q.t();
    
    // Rcpp::Rcout << "flag0 " << i << endl;
    
  }
  
  //Rcpp::Rcout << "flag2" << endl;
  
  arma::mat I = join_cols(join_rows(Ba,Bab),
                          join_rows(Bab.t(),Bb));
  I = I + S;
  
  arma::mat identmat = arma::eye(I.n_rows,I.n_cols);
  // arma::mat constmat = ones<mat>(I.n_rows,I.n_cols)*0.001;
  arma::mat Iinv = solve(I,identmat);
  arma::mat Irb = Iinv*S*Iinv.t();
  //
  arma::vec ASErb = sqrt(Irb.diag());
  arma::vec ASErb_alpha = ASErb.subvec(0,(num_class-1)*(x.n_cols+1)-1);
  arma::vec ASErb_beta = ASErb.subvec(ASErb.n_rows-num_feature*num_class*lbeta,ASErb.n_rows-1);
  //
  arma::vec ASEi = sqrt(Iinv.diag());
  arma::vec ASEi_alpha = ASEi.subvec(0,(num_class-1)*(x.n_cols+1)-1);
  arma::vec ASEi_beta = ASEi.subvec(ASEi.n_rows-num_feature*num_class*lbeta,ASEi.n_rows-1);
  
  Rcpp::List ret;
  // ret["I"] = I;
  // ret["S"] = S;
  // ret["q"] = qall;
  // ret["Irb"] = Irb;
  // ret["Ba"] = Ba;
  // ret["Bg"] = Bg;
  // ret["Bl"] = Bl;
  // ret["Bag"] = Bag;
  // ret["Bal"] = Bal;
  // ret["Bgl"] = Bgl;
  // ret["Blg"] = Blg;
  // ret["Bb"] = Bb;
  // ret["ASErb"] = ASErb;
  // ret["ASEi"] = ASEi;
  //ret["I"] = I;
  //ret["S"] = S;
  ret["alpharb"] = ASErb_alpha;
  ret["betarb"] = ASErb_beta;
  
  ret["alphai"] = ASEi_alpha;
  ret["betai"] = ASEi_beta;
  // ret["I"] = Irb;
  return ret;
}