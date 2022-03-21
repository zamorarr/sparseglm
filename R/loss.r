loss_exponential <- function(w, x, y) {
  z <- y * x
  loss_exponential2(w, z)
}

loss_exponential2 <- function(w, z) {
  mean(exp(-z %*% w))
}
