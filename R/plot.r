#' @export
plot.sparseglm_fit <- function(model, x, y) {
  p1 <- plot_cal(model, x, y)
  p2 <- plot_auc(model, x, y)
  p3 <- plot_acc(model, x, y)
  p4 <- plot_preds(model, x, y)
  gridExtra::grid.arrange(p1, p2, p3, p4)
}

plot_coefs <- function(model) {
  w <- coef(model, only_nz = TRUE, sort = TRUE)

  df <- tibble::tibble(w = w, name = names(w))
  df$name <- forcats::fct_reorder(df$name, df$w)

  ggplot2::ggplot(df, ggplot2::aes(x = name, y = w, fill = w >= 0)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal()
}

plot_cal <- function(model, x, y) {
  cal <- eval_cal(model, x, y, grouped = TRUE)
  cal_val <- eval_cal(model, x, y, grouped = FALSE)

  ggplot2::ggplot(cal, ggplot2::aes(x = predicted, y = observed)) +
    ggplot2::geom_abline(linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(size = n)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Risk Calibration",
      subtitle = sprintf("calibration error: %0.3f", cal_val),
      x = "predicted",
      y = "observed") +
    ggplot2::theme_minimal()
}

plot_auc <- function(model, x, y) {
  auc_val <- eval_auc(model, x, y)
  thresholds <- seq(0, 1, 0.1)
  df_prec_rec <- tibble::tibble(
    threshold = thresholds,
    precision = eval_precision(model, x, y, thresholds),
    recall = eval_recall(model, x, y, thresholds)
  )

  ggplot2::ggplot(df_prec_rec, ggplot2::aes(x = precision, y = recall, color = thresholds)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_viridis_c("threshold") +
    ggplot2::labs(title = "Precision vs Recall", subtitle = sprintf("auc: %0.3f", auc_val)) +
    ggplot2::theme_minimal()
}

plot_acc <- function(model, x, y) {
  # accuracy plot
  acc <- tibble::tibble(
    threshold = seq(0, 1, 0.1),
    acc = eval_acc(model, x, y, threshold = threshold))

  # calculate baseline accuracy
  acc_baseline <- mean(y == 1)
  if (acc_baseline < 0.5) acc_baseline <- 1 - acc_baseline

  ggplot2::ggplot(acc, ggplot2::aes(x = threshold, y = acc)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = acc_baseline, linetype = "dashed", color = "red") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous("accuracy", limits = c(0, 1)) +
    ggplot2::labs(
      title = "Accuracy",
      subtitle = sprintf("max: %0.02f (threshold = %0.02f)", max(acc$acc), acc$threshold[which.max(acc$acc)[1]])) +
    ggplot2::theme_minimal()
}

plot_preds <- function(model, x, y) {
  df <- tibble::tibble(y = y, prob = predict(model, x, type = "response"))
  ggplot2::ggplot(df, ggplot2::aes(x = prob)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = factor(y)),
      bins = 30,
      color = "white",
      position = "identity",
      #position = "dodge2",
      alpha = 0.5) +
    ggplot2::scale_fill_discrete("target", guide = "none") +
    ggplot2::labs(title = "Score Distribution") +
    ggplot2::facet_wrap(~y, ncol = 1) +
    ggplot2::theme_minimal()
}
