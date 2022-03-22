swap_features <- function(w, z, lambda) {
  cat("===swapping\n")
  # initialize support
  in_support <- w != 0

  # initialize failed swap counts
  swap_fails <- rep(0L, length(w))

  # initialize cost vector
  logh <- -z %*% w
  h <- exp(-z %*% w)
  cat(sprintf("  [initial] loss: %0.03f l0: %i\n", mean(h), sum(w != 0)))

  # loop until no more swaps
  while(TRUE) {
    # update support indices and sort by swap_fails
    support_idx <- which(in_support)
    support_idx <- support_idx[order(swap_fails[support_idx])]

    # loop over features in support index
    for (j in support_idx) {
      # try and swap
      out <- delete_or_swap(w, j, z, logh, lambda, in_support)

      if (out$is_swap) {
        # swap completed
        w <- out$w
        logh <- out$logh
        in_support <- out$in_support
        break
      } else {
        # swap failed
        swap_fails[j] <- swap_fails[j] + 1L
      }
    }

    # iterated over all features without a swap
    h <- exp(logh)
    H <- mean(h)
    cat(sprintf("  [final] loss: %0.03f l0: %i\n", H, sum(w != 0)))
    return(list(w = w, logh = logh))
  }
}

delete_or_swap <- function(w, j, z, logh, lambda, in_support) {
  # calculate current loss
  h <- exp(logh)
  H <- mean(h)
  l0 <- sum(in_support)
  wj <- w[j]

  # calculate auxillary variables
  neg_idx <- z[,j] == -1
  d <- sum(h[neg_idx])/sum(h)
  eps <- 1E-9
  d <- clamp(d, eps, 1 - eps)

  # try removing j
  w_zero <- w
  w_zero[j] <- 0
  H_zero <- H*(exp(wj)*(1 - d) + exp(-wj)*d)
  logh_zero <- logh
  logh_zero[neg_idx] <- logh_zero[neg_idx] + 0 - wj
  logh_zero[!neg_idx] <- logh_zero[!neg_idx] + wj - 0

  if (H_zero < (H + lambda)) {
    in_support[j] <- FALSE
    out <- list(w = w_zero, logh = logh_zero, in_support = in_support, is_swap = TRUE)
    return(out)
  }

  # try swapping j with something not in support?
  j_swaps <- which(!in_support)
  j_swaps <- sample(j_swaps, length(j_swaps))

  out_swaps <- lapply(j_swaps, function(j_swap) update_coef(w_zero[j_swap], j_swap, z[,j_swap], logh_zero, lambda))
  H_swaps <- vapply(out_swaps, function(out) mean(exp(out$logh)), double(1))
  best_idx <- which.min(H_swaps)[1]
  H_swap <- H_swaps[best_idx]

  if (H_swap < H) {
    wj_swap <- out_swaps[[best_idx]]$wj
    j_swap <- j_swaps[[best_idx]]

    w_swap <- w_zero
    w_swap[j_swap] <- wj_swap

    neg_idx <- z[,j_swap] == -1
    logh_swap <- logh_zero
    logh_swap[neg_idx] <- logh_swap[neg_idx] + wj_swap - 0
    logh_swap[!neg_idx] <- logh_swap[!neg_idx] + 0 - wj_swap

    in_support[j] <- FALSE
    in_support[j_swap] <- TRUE

    cat(sprintf("swapping %i with %i\n", j, j_swap))

    return(list(w = w_swap, logh = logh_swap, in_support = in_support, is_swap = TRUE))
  }

  # no removal, no swaps
  return(list(w = w, logh = logh, in_support = in_support, is_swap = FALSE))
}
