swap_features <- function(w, z) {
  # initialize support
  in_support <- w != 0

  # initialize failed swap counts
  swaps_fails <- rep(0L, length(w))

  # loop until no more swaps
  while(TRUE) {
    # update support indices and sort by swap_fails
    support_idx <- which(in_support)
    support_idx <- support_idx[order(swap_fails[support_idx])]

    # loop over features in support index
    for (j in support_idx) {
      # try and swap
      out <- delete_or_swap(w, j, z, in_support, lambda)

      if (out$is_swap) {
        # swap completed
        w <- out$w
        in_support <- out$in_support
        break
      } else {
        # swap failed
        swap_fails[j] <- swap_fails[j] + 1L
      }
    }

    # iterated over all features without a swap
    return(w)
  }
}

delete_or_swap <- function(w, j, z, in_support, lambda) {
  # calculate current loss
  loss_current <- loss_exponential2(w, z) + lambda*(sum(w != 0))

  # try removing j
  w_new <- w
  w_new[j] <- 0
  loss_remove <- loss_exponential2(w_new, z) + lambda*(sum(w_new != 0))

  if (loss_remove < loss_current) {
    w <- w_new
    in_support[j] <- FALSE
    out <- list(w = w, in_support = in_support, is_swap = TRUE)
    return(out)
  }

  # try swapping j with something not in support?

}

update_w <- function(x, y, s) {
  num_features <- ncol(x)
  model <- lm(y ~ . - 1, data = x[,s])
  w <- unname(coef(model))
  w2 <- double(num_features)
  w2[s] <- w
  return(w2)
}

delete_or_swap2 <- function(w, j, s, sc, x, y) {
  # calculate best current loss
  loss <- function(w, x, y) loss_logistic(w, x, y) + 0.001*sum(w^2)
  Lbest <- loss(w, x, y)

  # create new candidate solution by dropping feature j
  w2 <- w
  w2[j] <- 0

  # if loss is better w/o feature j
  if (loss(w2, x, y) <= Lbest) {
    # update support
    s2 <- s[s != j]

    # recalculate w2
    w2 <- update_w(x, y, s2)
    return(w2)
  }

  # if sc is empty, nothing else to check
  if (length(sc) < 1) return(w)

  # calcuate gradient on sc.
  grad <- grad_logistic(w, x, y)
  grad_sc <- grad[sc]
  grad_sc_order <- order(abs(grad_sc), decreasing = TRUE)

  for (j3 in grad_sc_order) {
    w3 <- linear_cut(w2, j3, Lbest)
    # if w and w2 are different
    if (any(w2 != w3)) {
      # update support
      s3 <- s[s != j]
      s3 <- sort(c(s, j3))

      # recalculate w
      w3 <- update_w(x, y, s3)
      return(w3)
    }
  }

  # return w
  w
}
