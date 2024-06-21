#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat LinProjCpp(arma::cube mu, arma::cube v, arma::mat gamma, 
                     arma::mat phi, arma::vec id, arma::vec obs, 
                     arma::vec time, arma::mat y, std::string cor){
  
  Rcpp::List lambda_res;
  Rcpp::List res;
  
  double num_class = mu.n_slices;
  double num_feature = mu.n_cols;
  arma::vec uniqueid = arma::unique(id);
  double n = uniqueid.n_elem;
  arma::mat LP = zeros<mat>(n,num_class);
  
  // Rcpp::Rcout << "span: " << span(num_feature) << std::endl;
  
  // Rcpp::Rcout << "iteration: " << i << " obsloglik: " << l << " diff: " << diff << endl;
  
  for (int i = 0; i < n ; ++i){
    
    uvec idx = find(id == i+1);
    // double numobs = arma::max(obs(idx));
    double numobs = idx.n_elem;
    
    arma::mat yi = y.rows(idx);
    arma::vec timei = time(idx);
    // Rcpp::Rcout << "flag 1 " << i+1 << "front " << idx.front() << "back " << idx.back() << "time " << timei << endl;
    
    // Rcpp::Rcout << "flag yi" << yi << endl;
    
    arma::cube vi = v.tube(span(idx.front(),idx.back()),span::all);
    arma::cube mui = mu.tube(span(idx.front(),idx.back()),span::all);
    arma::cube di = mui;
    
    if (numobs == 1){
      
      for (int c = 0; c < num_class ; ++c){
        
        di.slice(c) = yi - mui.slice(c);
        if (c > 0){
          LP(i,c) = accu((mui.slice(c) - mui.slice(0)) % (di.slice(c)/(vi.slice(c) % phi.row(c)) + di.slice(0)/(vi.slice(0) % phi.row(0))))/2;
        }
      }
      
    }else{
      
      arma::mat distime = zeros<mat>(numobs,numobs);
      // arma::vec distvec = arma::diff(timei);
      // arma::vec distvec = arma::ones(numobs-1);
      arma::vec distvec = arma::diff(obs(idx));
      
      for (int j = 0; j < numobs-1 ; ++j){
        // distime.submat(j,j+1,j,numobs-1) = arma::diff(timei,j+1);
        distime(j,span(j+1,numobs-1)) = cumsum(distvec.subvec(j,numobs-2)).t();
        distime(span(j+1,numobs-1),j) = cumsum(distvec.subvec(j,numobs-2));
      }
      // Rcpp::Rcout << "flag distime" << distime << endl;
      
      arma::sp_mat zcor0 = sp_mat(numobs*num_feature,numobs*num_feature);
      
      if (cor == "ar1"){
        for (int k = 0; k < num_feature; ++k){
          // Rcpp::Rcout << "flag distime" << pow(gamma(0,k),distime) << endl;
          for (int h1 = 0; h1 < numobs; ++h1){
            for (int h2 = 0; h2 < numobs; ++h2){
              zcor0(k*numobs+h1,k*numobs+h2) = std::pow(gamma(0,k),distime(h1,h2));
            }
          }
          // = distime.for_each([](mat::elem_type& val, mat::elem_type& gamma)  pow(val,gamma(0,k)); });    //arma::powmat(corrck,distime,1);
        }
      }else if (cor == "ind"){
        
        for (int k = 0; k < num_feature; ++k){
          zcor0 = eye(size(zcor0));
        }
        
      }
      
      // Rcpp::Rcout << "flag zcor0" << zcor0 << endl;
      
      for (int c = 0; c < num_class ; ++c){
        
        di.slice(c) = yi - mui.slice(c);
        
        // Rcpp::Rcout << "flag mui" << mui << endl;
        // Rcpp::Rcout << "flag vi" << vi << endl;
        if (c > 0){
          arma::sp_mat zcorc = sp_mat(numobs*num_feature,numobs*num_feature);
          
          if (cor == "ar1"){
            
            for (int k = 0; k < num_feature; ++k){
              // Rcpp::Rcout << "flag distime" << pow(gamma(0,k),distime) << endl;
              for (int h1 = 0; h1 < numobs; ++h1){
                for (int h2 = 0; h2 < numobs; ++h2){
                  zcorc(k*numobs+h1,k*numobs+h2) = std::pow(gamma(c,k),distime(h1,h2));
                }
              }
              // = distime.for_each([](mat::elem_type& val, mat::elem_type& gamma)  pow(val,gamma(0,k)); });    //arma::powmat(corrck,distime,1);
            }
          }
          else if (cor == "ind"){
            for (int k = 0; k < num_feature; ++k){
              zcorc = eye(size(zcorc));
            }
          }
          
          arma::mat mudvi0 = (mui.slice(c) - mui.slice(0))/sqrt(vi.slice(0).each_row() % phi.row(0));
          arma::mat mudvic = (mui.slice(c) - mui.slice(0))/sqrt(vi.slice(c).each_row() % phi.row(c));
          arma::mat vi0d0 = di.slice(0)/sqrt(vi.slice(0).each_row() % phi.row(0));
          arma::mat vicdc = di.slice(c)/sqrt(vi.slice(c).each_row() % phi.row(c));
          // Rcpp::Rcout << "flag mudvi0" << mudvi0 << endl;
          // Rcpp::Rcout << "flag vi0d0" << vi0d0 << endl;
          // Rcpp::Rcout << "flag mudvic" << mudvic << endl;
          // Rcpp::Rcout << "flag vicdc" << vicdc << endl;
          // Rcpp::Rcout << "flag mudvi0" << vectorise(mudvi0) << endl;
          
          LP(i,c) = accu( vectorise(mudvi0).t() * spsolve(zcor0,eye(size(zcor0)),"lapack") * vectorise(vi0d0) +
            vectorise(mudvic).t() * spsolve(zcorc,eye(size(zcorc)),"lapack") * vectorise(vicdc) )/2;
          
          // spsolve(zcor0,speye(size(zcor0)))
          
        }
      }
    }
  } 
  
  
  return exp(LP);
}