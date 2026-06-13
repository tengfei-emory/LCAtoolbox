#ifndef SLTCA_COMMON_H
#define SLTCA_COMMON_H

#include "RcppArmadillo.h"

// Pairwise wave distances |w_a - w_b| for one subject, used as the AR-1 exponent.
inline arma::mat sltca_distime(const arma::vec& waves){
  int m = waves.n_elem;
  arma::mat d(m, m);
  for (int a = 0; a < m; ++a){
    for (int b = 0; b < m; ++b){
      d(a, b) = std::fabs(waves(a) - waves(b));
    }
  }
  return d;
}

// Working correlation matrix for one feature, restricted to its observed rows
// (indices `keep` within the subject's rows). cor: "ind", "ar1" or "exch".
inline arma::mat sltca_cormat(const arma::mat& distime, const arma::uvec& keep,
                              double gamma, const std::string& cor){
  int m = keep.n_elem;
  arma::mat R = arma::eye(m, m);
  if (cor == "ar1"){
    for (int a = 0; a < m; ++a){
      for (int b = 0; b < m; ++b){
        if (a != b) R(a, b) = std::pow(gamma, distime(keep(a), keep(b)));
      }
    }
  }else if (cor == "exch"){
    for (int a = 0; a < m; ++a){
      for (int b = 0; b < m; ++b){
        if (a != b) R(a, b) = gamma;
      }
    }
  }
  return R;
}

// GEE quasi-score (qseg, length lbeta) and information block (Dblock,
// lbeta x lbeta) of one feature for one subject and class. Rows where the
// feature is missing (non-finite y) are dropped, so each feature contributes
// through its own observed visits only.
inline void sltca_qscore(const arma::vec& yjl, const arma::vec& mujlc,
                         const arma::vec& t0j, const arma::vec& geex,
                         const std::string& dist, double phi_lc, double gamma_lc,
                         const std::string& cor, const arma::mat& distime,
                         int lbeta, arma::vec& qseg, arma::mat& Dblock){

  qseg.zeros(lbeta);
  Dblock.zeros(lbeta, lbeta);

  arma::uvec keep = arma::find_finite(yjl);
  int m = keep.n_elem;
  if (m == 0) return;

  arma::vec yo = yjl(keep);
  arma::vec muo = mujlc(keep);
  arma::vec to = t0j(keep);

  // vlink = sqrt(V(mu)); dscale = dmu/deta under the canonical link
  arma::vec vlink(m), dscale(m);
  if (dist == "normal"){
    vlink.ones();
    dscale.ones();
  }else if (dist == "poi"){
    vlink = arma::sqrt(muo);
    dscale = muo;
  }else{ // "bin"
    vlink = arma::sqrt(muo % (1 - muo));
    dscale = muo % (1 - muo);
  }

  arma::mat dmu(m, lbeta);
  dmu.col(0) = dscale;
  dmu.col(1) = dscale % to;
  for (int g = 2; g < lbeta; ++g){
    dmu.col(g) = dscale * geex(g - 2);
  }

  arma::mat R = sltca_cormat(distime, keep, gamma_lc, cor);
  arma::mat V = phi_lc * ((vlink * vlink.t()) % R); // diag(vlink) R diag(vlink)

  qseg = dmu.t() * arma::solve(V, yo - muo);
  Dblock = dmu.t() * arma::solve(V, dmu);
}

#endif
