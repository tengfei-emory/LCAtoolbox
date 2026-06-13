#include "RcppArmadillo.h"
#include "sltca_common.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// Linear projection of the approximated likelihood ratio (class c vs class 1).
// The working covariance is block diagonal across features (conditional
// independence between features given the class), so the quadratic form is
// accumulated feature by feature; for each feature only its observed
// (finite) y entries contribute.
// [[Rcpp::export]]
arma::mat LinProjCpp(arma::cube mu, arma::cube v, arma::mat gamma,
                     arma::mat phi, arma::vec id, arma::vec obs,
                     arma::vec time, arma::mat y, std::string cor){

  int num_class = mu.n_slices;
  int num_feature = mu.n_cols;
  arma::vec uniqueid = arma::unique(id);
  int n = uniqueid.n_elem;
  arma::mat LP = zeros<mat>(n, num_class);

  for (int i = 0; i < n; ++i){

    uvec idx = find(id == i + 1);
    arma::mat yi = y.rows(idx);
    arma::mat distime = sltca_distime(obs(idx));

    arma::mat mu0i = mu.slice(0).rows(idx);
    arma::mat v0i = v.slice(0).rows(idx);

    for (int c = 1; c < num_class; ++c){

      arma::mat muci = mu.slice(c).rows(idx);
      arma::mat vci = v.slice(c).rows(idx);

      double lp = 0;

      for (int k = 0; k < num_feature; ++k){

        arma::vec yk = yi.col(k);
        arma::uvec keep = find_finite(yk);
        if (keep.n_elem == 0) continue;

        arma::vec yo = yk(keep);
        arma::vec m0 = mu0i.col(k); m0 = m0(keep);
        arma::vec mc = muci.col(k); mc = mc(keep);
        arma::vec s0 = sqrt(v0i.col(k) * phi(0, k)); s0 = s0(keep);
        arma::vec sc = sqrt(vci.col(k) * phi(c, k)); sc = sc(keep);

        arma::mat R0 = sltca_cormat(distime, keep, gamma(0, k), cor);
        arma::mat Rc = sltca_cormat(distime, keep, gamma(c, k), cor);

        lp += dot((mc - m0) / s0, solve(R0, (yo - m0) / s0)) +
              dot((mc - m0) / sc, solve(Rc, (yo - mc) / sc));
      }

      LP(i, c) = lp / 2;
    }
  }

  return exp(LP);
}
