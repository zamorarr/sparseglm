coord_descent <- function(w, z, lambda, max_iters = 100) {
  cat("===coord descent\n")

  # negative values of zj (the margin)
  neg_idx <- z == -1

  # initialize cost vector
  logh <- -z %*% w

  # calculate initial loss
  h <- exp(logh)
  H <- mean(h)
  l0 <- sum(w != 0)
  loss <- H + lambda*l0
  cat(sprintf("  [initial] loss: %0.03f l0: %i\n", loss, l0))

  # update coefficients
  for (iter in seq_len(max_iters)) {
    loss_old <- loss

    # loop over every feature
    for (j in seq_along(w)) {
      delta_j <- update_coef(w[j], j, neg_idx[,j], logh, lambda)
      w[j] <- w[j] + delta_j
      logh <- update_logh(logh, neg_idx[,j], delta_j)
    }

    # calculate loss
    H <- mean(exp(logh))
    l0 <- sum(w != 0)
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
  list(w = w, logh = logh)
}

update_coef <- function(wj, j, neg_idx, logh, lambda) {
  # current total cost
  h <- exp(logh)
  H <- mean(h)

  # calculate auxillary variables
  d <- sum(h[neg_idx])/sum(h)
  eps <- 1E-9
  d <- clamp(d, eps, 1 - eps)

  # candidate update
  wj_new <- wj + 0.5*(log(1 - d) - log(d))
  H_new <- 2*H*sqrt(d*(1 - d))

  # compare to setting wj -> 0
  wj_zero <- 0
  H_zero <- H*(exp(wj)*(1 - d) + exp(-wj)*d)

  # check if it is good enough to update
  loss_new <- H_new + lambda*(wj_new != 0)
  loss_zero <- H_zero

  if (loss_new < loss_zero) {
    wj_best <- wj_new
  } else {
    wj_best <- wj_zero
  }

  # return change
  wj_best - wj
}

update_logh <- function(logh, neg_idx, delta) {
  logh[neg_idx] <- logh[neg_idx] + delta
  logh[!neg_idx] <- logh[!neg_idx] - delta
  logh
}
