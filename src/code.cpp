#include <RcppArmadillo.h>
#include "cdexponential.h"
//using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List coord_descent_cpp(arma::vec &w, const arma::mat &z, double lambda, uint max_iters) {
  CDExponential cd(w, z, lambda);
  return cd.fit(max_iters);
}
