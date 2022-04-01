#include "cdexponential.h"

uint n_nzcoef(arma::vec &w) {
  // don't count intercept
  return arma::accu(w != 0) - (w(0) != 0);
}

double clamp(double x, double low, double high) {
  return std::min(std::max(x, low), high);
}

double compute_d(const arma::vec &h, arma::uvec indices) {
  double d = arma::sum(h.elem(indices))/arma::sum(h);
  return clamp(d, 1E-9, 1 - 1E-9);
}

// constructor
CDExponential::CDExponential(arma::vec &_w, const arma::mat &_z, double _lambda):
  w(_w), z(_z), lambda(_lambda) {

  // initialize
  for (uint j = 0; j < z.n_cols; j++) {
    indices[j] = arma::find(z.col(j) == -1);
  }

  h = arma::exp(-z * w);
  H = arma::mean(h);
}

Rcpp::List CDExponential::fit(uint max_iters) {
  uint l0 = n_nzcoef(w);
  double loss =  H + lambda*l0;
  Rcpp::Rcout << "  [initial] loss: " << loss << " l0: " << l0 << std::endl;

  // loop max_iter times
  double loss_old;
  double abs_tol;
  double rel_tol;
  bool is_done = false;

  // iterate max_iters times at most
  for (uint iter = 0; iter < max_iters; iter++) {
    loss_old = loss;

    // iterate over every feature
    for (uint j = 0; j < w.size(); j++) {
      update_coef(j);
    }

    // compute new loss
    l0 = n_nzcoef(w);
    loss = H + lambda*l0;

    // check if we can break
    abs_tol = loss_old - loss;
    rel_tol = 1 - loss/loss_old;
    is_done = abs_tol <= 1E-9 || rel_tol <= 1E-6;
    if (is_done) break;
  }


  Rcpp::Rcout << "  [final] loss: " << loss << " l0: " << l0 << std::endl;

  // tidy results
  Rcpp::NumericVector w_vec(w.size());
  for (uint i = 0; i < w.size(); i++) {
    w_vec(i) = w(i);
  }

  Rcpp::NumericVector h_vec(h.size());
  for (uint i = 0; i < h.size(); i++) {
    h_vec(i) = h(i);
  }

  return Rcpp::List::create(
    Rcpp::Named("w") = w_vec,
    Rcpp::Named("h") = h_vec);
}

void CDExponential::update_h(uint j, double delta) {
  h %= arma::exp(-delta*z.col(j));
}

double CDExponential::compute_H(double H, double d, double delta) {
  return H*(std::exp(-delta)*(1 - d) + std::exp(delta)*d);
}

void CDExponential::update_coef(uint j) {
  // calculate auxillary variables
  double d = compute_d(h, indices[j]);

  // candidate update
  double delta_new = 0.5*(log(1 - d) - log(d));
  double H_new = 2*H*sqrt(d*(1 - d));

  // compare to setting wj -> 0
  double delta_zero = -w(j);
  double H_zero = compute_H(H, d, delta_zero);

  // check if it is good enough to update
  double loss_new = H_new + lambda*(j != 0);
  double loss_zero = H_zero;

  // return best
  if (loss_new < loss_zero) {
    update_h(j, delta_new);
    H = H_new;
    w(j) += delta_new;
  } else {
    update_h(j, delta_zero);
    H = H_zero;
    w(j) += delta_zero;
  }
}

