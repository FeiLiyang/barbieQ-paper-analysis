library(knitr)

# Specific HTML settings for files under subdirectory of ./analysis
html_settings_subdirectory <- list(
  toc = TRUE,
  toc_depth = 4,
  toc_float = list(collapsed = TRUE, smooth_scroll = TRUE),
  number_sections = TRUE,
  theme = "cosmo",
  highlight = "textmate"
)

# Determine the subdirectory from the Rmd path
set_html_settings <- function(rmd_path) {
  if (grepl("analysis/WuC", rmd_path)) {
    opts_knit$set(output_dir = "../../docs/WuC", html_settings = html_settings_subdirectory)
  } else {
    opts_knit$set(output_dir = "../../docs", html_settings = NULL)
  }
}

# Example usage (adapt this in your Rmd files):
# set_html_settings(knitr::current_input())
