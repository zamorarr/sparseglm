coord_descent <- function(z, lambda, max_iters = 100, w_init = NULL) {
  n <- nrow(z)
  p <- ncol(z)

  # initialize coefficient vector
  if (is.null(w_init)) {
    w <- runif(p, min = -1, max = 1)
  } else {
    w <- w_init
  }

  # initialize cost vector
  h <- exp(-z %*% w)
  cat(sprintf("[initial] loss: %0.03f l0: %i\n", mean(h), sum(w != 0)))

  # update coefficients
  for (iter in seq_len(max_iters)) {
    w_old <- w

    # loop over every feature
    for (j in seq_along(w)) {
      out <- update_coef(w[j], j, z[,j], h, lambda)
      w[j] <- out$wj
      h <- out$h
    }

    # check if we can break
    max_diff <- max(abs(w - w_old))
    if (max_diff <= 1E-8) break
  }

  cat(sprintf("[final] %0.03f loss: %0.03f l0: %i\n", max_diff, mean(h), sum(w != 0)))

  # return results
  w
}

update_coef <- function(wj, j, zj, h, lambda) {
  #zj <- z[,j]
  neg_idx <- zj == -1

  # current total cost
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
  H_zero <- H + (1 - h[j])/length(h)
  #H_zero <- mean(h[-j]) # this is not correct

  # check if it is good enough to update
  #loss <- H + lambda*(wj != 0)
  loss_new <- H_new + lambda*(wj_new != 0)
  loss_zero <- H_zero

  if (loss_new < loss_zero) {
    wj_best <- wj_new
  } else {
    wj_best <- wj_zero
  }

  # update h
  h[neg_idx] <- h[neg_idx]*exp(wj_best - wj)
  h[!neg_idx] <- h[!neg_idx]*exp(wj - wj_best)

  # return results
  list(wj = wj_best, h = h)
}
