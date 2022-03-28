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
