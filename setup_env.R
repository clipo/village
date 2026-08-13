#' Environment setup for the local Susquehanna contemporaneity analysis.
#'
#' Installs the R packages required by the analysis into the project renv
#' library, then installs the CmdStan toolchain. Safe to re-run; already
#' installed components are skipped.

options(repos = c(
  stan = "https://stan-dev.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

cran_pkgs <- c(
  "yaml", "readxl", "dplyr", "tidyr", "purrr", "stringr",
  "ggplot2", "scales", "patchwork", "ggridges",
  "systemfonts", "ragg",
  "rcarbon", "posterior", "loo", "bayesplot",
  "cmdstanr"
)

missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  cat("Installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, quiet = FALSE)
} else {
  cat("All R packages already present.\n")
}

failed <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed)) {
  cat("\nFAILED TO INSTALL:", paste(failed, collapse = ", "), "\n")
} else {
  cat("\nAll R packages installed.\n")
}

# CmdStan toolchain -----------------------------------------------------------
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver)) {
    cat("\nInstalling CmdStan (this takes several minutes)...\n")
    cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
    cmdstanr::install_cmdstan(cores = max(1L, parallel::detectCores() - 2L),
                              overwrite = FALSE)
    ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  }
  cat("\nCmdStan version:", if (is.null(ver)) "NOT AVAILABLE" else ver, "\n")
}

cat("\nSETUP COMPLETE\n")
