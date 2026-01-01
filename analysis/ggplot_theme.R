# Set ggplot theme

## ggplot theme
strip_zeros <- function(x) {
  x <- format(x, scientific = FALSE)
  x <- sub("0+$", "", x)
  sub("\\.$", "", x)
}
## ggplot theme
theme_barbieQ <- theme_barbie <- function(base_size = 12, base_family = "") {
  theme_linedraw(base_size = base_size, base_family = base_family) +
    theme(
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = base_size * 1.1),
      axis.text  = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size * 0.9),
      legend.key = element_blank(),
      plot.title = element_text(size = base_size * 1.1),
      plot.subtitle = element_text(size = base_size)
    )
}

theme_barbieQ_size <- theme_barbie <- function(base_size = 12) {
  theme(
    axis.title = element_text(size = base_size * 1.1),
    axis.text  = element_text(size = base_size),
    legend.title = element_text(size = base_size),
    legend.text  = element_text(size = base_size * 0.9),
    plot.title = element_text(size = base_size * 1.1),
    plot.subtitle = element_text(size = base_size)
  )
}
