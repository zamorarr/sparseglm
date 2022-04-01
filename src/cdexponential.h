#ifndef CDEXPONENTIAL_H
#define CDEXPONENTIAL_H

#include <RcppArmadillo.h>

class CDExponential {
private:
  // provided
  arma::vec w; // coefficients to optimize {p x 1}
  arma::mat z; // data {n x p}
  double lambda; // regularization constant

  // derived
  std::unordered_map<uint, arma::uvec> indices; // indices where z == -1
  arma::vec h; // {n x 1}
  double H;

public:
  CDExponential(arma::vec &_w, const arma::mat &_z, double _lambda);
  Rcpp::List fit(uint max_iters);
  void update_coef(uint j);
  void update_h(uint j, double delta);
  double compute_H(double H, double d, double delta);
};

#endif
