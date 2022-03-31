# my swap loss returned is not matching up with mean(exp(-z %*% w))?
swap_features <- function(w, z, lambda) {
  cat("===swapping\n")
  # initialize support
  in_support <- w != 0

  # initialize failed swap counts
  swap_fails <- rep(0L, length(w))

  # initialize cost vector
  logh <- -z %*% w
  h <- exp(logh)

  # calculate loss
  H <- mean(h)
  l0 <- n_nzcoef(w)
  loss <- H + lambda*l0
  cat(sprintf("  [initial] loss: %0.03f l0: %i\n", loss, l0))

  # loop until no more swaps
  while(TRUE) {
    # update support indices and sort by swap_fails
    support_idx <- which(in_support)
    support_idx <- support_idx[order(swap_fails[support_idx])]

    # loop over features in support index
    for (j in support_idx) {
      # never swap the intercept
      if (j == 1) next

      # try and swap
      out <- delete_or_swap(w, j, z, h, H, lambda, in_support)

      if (out$is_swap) {
        # swap completed
        w <- out$w
        h <- out$h
        H <- out$H
        in_support <- which(w != 0)
        break
      } else {
        # swap failed
        swap_fails[j] <- swap_fails[j] + 1L
      }
    }

    # iterated over all features without a swap
    l0 <- n_nzcoef(w)
    loss <- H + lambda*l0
    cat(sprintf("  [final] loss: %0.03f l0: %i\n", loss, l0))
    return(list(w = w, h = h, H = H))
  }
}

delete_or_swap <- function(w, j, z, h, H, lambda, in_support) {
  # calculate current loss
  l0 <- n_nzcoef(w)
  wj <- w[j]

  # calculate auxillary variables
  neg_idx <- z[,j] == -1
  d <- compute_d(h, neg_idx)

  # try removing j
  w_zero <- w
  w_zero[j] <- 0
  H_zero <- update_H(H, d, -wj)
  h_zero <- update_h(h, neg_idx, -wj)

  if (H_zero < (H + lambda)) {
    cat(sprintf("  zeroing %i\n", j))
    out <- list(w = w_zero, h = h_zero, H = H_zero, is_swap = TRUE)
    return(out)
  }

  # try swapping j with something not in support?
  j_swaps <- which(!in_support)
  out_swaps <- lapply(j_swaps, function(j) {
    wj <- w_zero[j]
    neg_idx <- z[,j] == -1

    # update coef
    update_coef(wj, neg_idx, h_zero, H_zero, lambda)
  })

  # find best swap
  H_swaps <- vapply(out_swaps, function(out) out$H, double(1))
  best_idx <- which.min(H_swaps)[1]
  H_swap <- H_swaps[best_idx]

  # swap has lower loss
  if (H_swap < H) {
    j_swap <- j_swaps[[best_idx]]
    w_swap <- w_zero
    w_swap[j_swap] <- out_swaps[[best_idx]]$wj
    h_swap <- out_swaps[[best_idx]]$h

    cat(sprintf("  swapping %i with %i\n", j, j_swap))
    return(list(w = w_swap, h = h_swap, H = H_swap, is_swap = TRUE))
  }

  # no removal, no swaps
  return(list(w = w, h = h, H = H, is_swap = FALSE))
}
