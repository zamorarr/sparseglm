#' @importFrom stats predict
#' @export
stats::predict

#' @export
predict.sparseglm_fit <- function(model, x, type = c("link", "response")) {
  type <- match.arg(type)

  score <- x %*% model$w
  score <- score[,1]
  if (identical(type, "response")) score_to_prob(score)
  else score
}

score_to_prob <- function(score) {
  1/(1 + exp(-2*score))
}
