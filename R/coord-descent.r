coord_descent <- function(z, lambda, max_iters = 100, w_init = NULL) {
  cat("===coord descent\n")
  n <- nrow(z)
  p <- ncol(z)

  # initialize coefficient vector
  if (is.null(w_init)) {
    w <- runif(p, min = -1, max = 1)
  } else {
    w <- w_init
  }

  # initialize cost vector
  logh <- -z %*% w
  h <- exp(-z %*% w)
  cat(sprintf("  [initial] loss: %0.03f l0: %i\n", mean(h), sum(w != 0)))

  # update coefficients
  for (iter in seq_len(max_iters)) {
    w_old <- w

    # loop over every feature
    for (j in seq_along(w)) {
      out <- update_coef(w[j], j, z[,j], logh, lambda)
      w[j] <- out$wj
      logh <- out$logh
    }

    # check if we can break
    max_diff <- max(abs(w - w_old))
    if (max_diff <= 1E-8) break
  }

  h <- exp(logh)
  H <- mean(h)
  cat(sprintf("  [final] loss: %0.03f l0: %i\n", H, sum(w != 0)))

  # return results
  list(w = w, logh = logh)
}

update_coef <- function(wj, j, zj, logh, lambda) {
  neg_idx <- zj == -1

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

  # update h
  logh[neg_idx] <- logh[neg_idx] + wj_best - wj
  logh[!neg_idx] <- logh[!neg_idx] + wj - wj_best

  # return results
  list(wj = wj_best, logh = logh)
}
