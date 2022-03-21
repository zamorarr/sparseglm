#' Create SparseGLM model
sparseglm <- function(x, y, lambda = 1E-3, max_iters = 100) {
  # initialize
  w_swap <- NULL
  z <- y*x

  # 1. coordinate descent
  w <- coord_descent(z, lambda, max_iters = max_iters, w_init = w_swap)
  in_sup <- w != 0

  # 2. swap features
  #w_swap <- swap_features(w, z, lambda)

  # return output object
  output <- list(w = w, in_sup = in_sup, names = colnames(x), lambda = lambda)
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
