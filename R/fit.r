#' Create SparseGLM model
sparseglm <- function(x, y, lambda = 1E-3, max_iters = 100) {
  # initialize
  w_swap <- NULL
  z <- y*x

  while(TRUE) {
    # 1. coordinate descent
    out <- coord_descent(z, lambda, max_iters = max_iters, w_init = w_swap)
    w <- out$w
    logh <- out$logh

    # 2. swap features
    out <- swap_features(w, z, lambda)
    w_swap <- out$w
    logh_swap <- out$logh

    # 3. check stopping criteria
    loss <- mean(exp(logh)) + lambda*(sum(w != 0))
    loss_swap <- mean(exp(logh_swap)) + lambda*(sum(w_swap != 0))

    #max_diff <- max(abs(w - w_swap))
    #is_done <- max_diff <= 1E-7
    is_done <- loss - loss_swap <= 1E-3
    if (is_done) break
    #print(w - w_swap)
  }

  # return output object
  output <- list(w = w, names = colnames(x), in_support = w != 0, l0 = sum(w != 0), lambda = lambda)
  structure(output, class = c("sparseglm_fit", class(x)))
}

#' @export
coef.sparseglm_fit <- function(object, only_nz = FALSE, sort = FALSE) {
  w <- setNames(object$w, object$names)

  if (only_nz) {
    w <- w[object$in_sup]
  }

  if (sort) {
    w <- w[order(w, decreasing = TRUE)]
  }

  w
}
