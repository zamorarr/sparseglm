eval_auc <- function(model, x, y) {
  ind_pos <- which(y > 0)
  ind_neg <- which(y <= 0)

  n_pos <- length(ind_pos)
  n_neg <- length(ind_neg)

  s <- predict(model, x, type = "link")
  auc <- vapply(ind_pos, function(i) sum(s[i] > s[ind_neg]), integer(1))
  auc <- sum(auc)

  exp(log(auc) - log(n_pos) - log(n_neg))
}

eval_acc <- function(model, x, y, threshold = 0.5) {
  probs <- predict(model, x, type = "response")

  vapply(threshold, function(r) {
    y_est <- 2L*(probs >= r) - 1L
    mean(y_est == y)
  }, double(1))
}

eval_precision <- function(model, x, y, threshold = 0.5) {
  probs <- predict(model, x, type = "response")

  vapply(threshold, function(r) {
    y_est <- 2L*(probs >= r) - 1L
    idx <- which(y_est == 1)
    n <- length(idx)

    if (n == 0) return(1)

    n_correct <- sum(y[idx] == 1)
    n_correct/n
  }, double(1))
}

eval_recall <- function(model, x, y, threshold = 0.5) {
  probs <- predict(model, x, type = "response")

  vapply(threshold, function(r) {
    y_est <- 2L*(probs >= r) - 1L
    idx <- which(y == 1)
    n <- length(idx)

    n_correct <- sum(y_est[idx] == 1)
    n_correct/n
  }, double(1))
}
