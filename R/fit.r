#' Create SparseGLM model
sparseglm <- function(x, y, lambdas = 10^seq(-3, -1), max_iters = 100) {
  z <- y*x
  lambdas <- sort(lambdas, decreasing = TRUE)

  out <- vector("list",length = length(lambdas))
  w_init <- rep(0, ncol(x))

  for (i in seq_along(lambdas)) {
    cat("lambda: ", lambdas[i], "\n")
    out[[i]] <- sparseglm_fit(z, lambda = lambdas[i], max_iters = max_iters, w_init = w_init)
    w_init <- out[[i]]$w
  }

  out
}

sparseglm_fit <- function(z, lambda, max_iters = 100, w_init = NULL) {
  # initialize coefficient vector
  if (is.null(w_init)) {
    w_swap <- runif(ncol(z), min = -1, max = 1)
  } else {
    w_swap <- w_init
  }

  while(TRUE) {
    # 1. coordinate descent
    out <- coord_descent(w_swap, z, lambda, max_iters = max_iters)
    w <- out$w
    logh <- out$logh
    loss <- mean(exp(logh)) + lambda*(sum(w != 0))

    # 2. swap features
    out <- swap_features(w, z, lambda)
    w_swap <- out$w
    logh_swap <- out$logh
    loss_swap <- mean(exp(logh_swap)) + lambda*(sum(w_swap != 0))

    # 3. check stopping criteria
    abs_tol <- loss - loss_swap
    rel_tol <- 1 - loss_swap/loss
    cat(sprintf("  reltol: %0.6f abstol: %0.9f\n", rel_tol, abs_tol))
    is_done <- abs_tol <= 1E-9 || rel_tol <= 1E-6
    if (is_done) break
  }

  # return output object
  output <- list(w = w, names = colnames(z), loss = loss, l0 = sum(w != 0), lambda = lambda)
  structure(output, class = c("sparseglm_fit", class(output)))
}

#' @export
coef.sparseglm_fit <- function(object, only_nz = FALSE, sort = FALSE) {
  w <- setNames(object$w, object$names)

  if (only_nz) {
    w <- w[w != 0]
  }

  if (sort) {
    w <- w[order(w, decreasing = TRUE)]
  }

  w
}

summarize_models <- function(models, x, y) {
  df <- tibble::as_tibble(purrr::transpose(models))
  df <- tidyr::unnest(df, c(loss, l0, lambda))
  df$auc <- vapply(models, function(m) eval_auc(m, x, y), double(1))
  df$acc <- vapply(models, function(m) eval_acc(m, x, y), double(1))
  df$prec <- vapply(models, function(m) eval_precision(m, x, y), double(1))
  df$recall <- vapply(models, function(m) eval_recall(m, x, y), double(1))
  df
}
