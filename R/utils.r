clamp <- function(x, a, b) {
  x[x < a] <- a
  x[x > b] <- b
  x
}
