coord_descent <- function(w, z, lambda, max_iters = 100) {
  cat("===coord descent\n")

  # negative values of zj (the margin)
  neg_idx <- z == -1

  # initialize cost vector
  logh <- -z %*% w
  h <- exp(logh)

  # calculate initial loss
  H <- mean(h)
  l0 <- n_nzcoef(w)
  loss <- H + lambda*l0
  cat(sprintf("  [initial] loss: %0.03f l0: %i\n", loss, l0))

  # update coefficients
  for (iter in seq_len(max_iters)) {
    loss_old <- loss

    # loop over every feature
    for (j in seq_along(w)) {
      wj <- w[j]
      neg_j <- neg_idx[,j]
      out <- update_coef(wj, neg_j, h, H, lambda, is_intercept = j == 1)

      w[j] <- out$wj
      h <- out$h
      H <- out$H
    }

    # calculate loss
    l0 <- n_nzcoef(w)
    loss <- H + lambda*l0

    # check if we can break
    abs_tol <- loss_old - loss
    rel_tol <- 1 - loss/loss_old
    is_done <- abs_tol <= 1E-9 || rel_tol <= 1E-6
    if (is_done) break
  }

  cat(sprintf("  [final] loss: %0.03f l0: %i\n", loss, l0))
  cat(sprintf("  finished at %i iterations\n", iter))

  # return results
  list(w = w, h = h)
}

update_coef <- function(wj, neg_idx, h, H, lambda, is_intercept = FALSE) {
  # current total cost
  #H <- mean(h)

  # calculate auxillary variables
  d <- compute_d(h, neg_idx)

  # candidate update
  delta_new <- 0.5*(log(1 - d) - log(d))
  H_new <- 2*H*sqrt(d*(1 - d))

  # compare to setting wj -> 0
  delta_zero <- -wj
  H_zero <- update_H(H, d, delta_zero)

  # check if it is good enough to update
  loss_new <- H_new + lambda*(!is_intercept)
  loss_zero <- H_zero

  # return best
  if (loss_new < loss_zero) {
    h_new <- update_h(h, neg_idx, delta_new)
    list(wj = wj + delta_new, H = H_new, h = h_new)
  } else {
    h_zero <- update_h(h, neg_idx, delta_zero)
    list(wj = wj + delta_zero, H = H_zero, h = h_zero)
  }
}

compute_d <- function(h, neg_idx) {
  #d <- sum(h[neg_idx])/sum(h)
  d <- sum(h*neg_idx)/sum(h)
  eps <- 1E-9
  d <- clamp(d, eps, 1 - eps)
}

update_H <- function(H, d, delta) {
  H*(exp(-delta)*(1 - d) + exp(delta)*d)
}

update_logh <- function(logh, neg_idx, delta) {
  #logh[neg_idx] <- logh[neg_idx] + delta
  #logh[!neg_idx] <- logh[!neg_idx] - delta
  logh + delta*(2*neg_idx - 1)
}

update_h <- function(h, neg_idx, delta) {
  h*exp(delta*(2*neg_idx - 1))
}

n_nzcoef <- function(w) sum(w[-1] != 0)

