# packages.R
# Install all dependencies for LipiDash Shiny app
# R 4.5.2 | Bioconductor 3.22
# Run during Docker image build; do NOT source at runtime.

# ── Helpers ────────────────────────────────────────────────────────────────────
install_if_missing <- function(pkgs, installer = install.packages, ...) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    installer(to_install, ...)
  } else {
    message("All packages already installed: ", paste(pkgs, collapse = ", "))
  }
}

# ── CRAN packages ──────────────────────────────────────────────────────────────
cran_packages <- c(
  # Shiny & UI
  "shiny",
  "DT",
  "shinyWidgets",
  "bslib",
  "htmltools",
  "htmlwidgets",
  "httpuv",
  "fontawesome",
  "crosstalk",
  "promises",
  "later",
  "jquerylib",
  "sass",
  "sourcetools",
  "xtable",

  # Data wrangling (tidyverse core)
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "purrr",
  "stringr",
  "rlang",
  "generics",
  "magrittr",
  "glue",
  "vctrs",
  "pillar",
  "pkgconfig",
  "lifecycle",
  "tidyselect",
  "vroom",
  "hms",
  "tzdb",
  "bit",
  "bit64",
  "cli",
  "utf8",
  "backports",
  "withr",
  "fs",

  # Plotting
  "ggplot2",
  "plotly",
  "ggpubr",
  "ggrepel",
  "ggsci",
  "ggsignif",
  "factoextra",
  "pheatmap",
  "svglite",
  "ragg",
  "viridis",
  "viridisLite",
  "RColorBrewer",
  "colorspace",
  "scales",
  "gtable",
  "isoband",
  "farver",
  "labeling",
  "gridExtra",
  "cowplot",
  "systemfonts",
  "textshaping",

  # Statistics & multivariate
  "FactoMineR",
  "rstatix",
  "car",
  "carData",
  "emmeans",
  "estimability",
  "multcompView",
  "lme4",
  "lmtest",
  "pbkrtest",
  "nloptr",
  "minqa",
  "nlme",
  "mgcv",
  "Matrix",
  "MatrixModels",
  "SparseM",
  "quantreg",
  "MASS",
  "boot",
  "survival",
  "broom",
  "doBy",
  "reformulas",
  "mvtnorm",
  "leaps",
  "scatterplot3d",
  "flashClust",
  "cluster",
  "ellipse",
  "corrplot",
  "dendextend",
  "polynom",

  # Network
  "visNetwork",

  # Time series / forecasting (indirect deps)
  "forecast",
  "timeDate",
  "zoo",
  "fracdiff",
  "urca",
  "tseries",

  # Misc / infrastructure
  "R6",
  "Rcpp",
  "RcppArmadillo",
  "RcppEigen",
  "RcppTOML",
  "BH",
  "cpp11",
  "abind",
  "modelr",
  "data.table",
  "fastmatch",
  "fastmap",
  "memoise",
  "cachem",
  "digest",
  "jsonlite",
  "yaml",
  "knitr",
  "rmarkdown",
  "xfun",
  "evaluate",
  "highr",
  "formatR",
  "tinytex",
  "progress",
  "prettyunits",
  "rbibutils",
  "Rdpack",
  "openssl",
  "curl",
  "httr",
  "jose",
  "sys",
  "askpass",
  "PKI",
  "base64enc",
  "commonmark",
  "mime",
  "codetools",
  "lattice",
  "lazyeval",
  "numDeriv",
  "microbenchmark",
  "nnet",
  "rappdirs",
  "rstudioapi",
  "clipr",
  "otel",
  "packrat",
  "renv",
  "rsconnect",
  "snowflakeauth",
  "snow",
  "lambda.r",
  "futile.logger",
  "futile.options",
  "S7",
  "Formula",
  "Deriv",
  "sbw",          # dependency of some packages
  "BiocManager"
)

install_if_missing(
  cran_packages,
  installer = install.packages,
  repos     = "https://cran.rstudio.com",
  Ncpus     = parallel::detectCores(),
  quiet     = TRUE
)

# ── Bioconductor packages ──────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com")
}

BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

bioc_packages <- c(
  "BiocParallel",   # parallel backend used by fgsea
  "BiocVersion",    # locks Bioc version
  "fgsea"           # fast gene set enrichment (used in Enrichment tab)
)

install_if_missing(
  bioc_packages,
  installer = function(pkgs, ...) BiocManager::install(pkgs, ask = FALSE, update = FALSE, ...),
  Ncpus     = parallel::detectCores()
)

message("All packages installed successfully.")
