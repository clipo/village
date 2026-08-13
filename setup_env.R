#' Environment setup for the local Susquehanna contemporaneity analysis.
#'
#' Installs the R packages required by the analysis into the project renv
#' library, then installs the CmdStan toolchain. Safe to re-run; already
#' installed components are skipped.

options(repos = c(
  stan = "https://stan-dev.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# sf (an rcarbon dependency) resolves gdal-config from PATH. Anaconda ships a
# gdal-config that does not match the Homebrew GDAL libraries, so put Homebrew
# ahead of it. Without this, sf fails to build and rcarbon fails with it.
Sys.setenv(PATH = paste(
  "/opt/homebrew/opt/gdal/bin",
  "/opt/homebrew/opt/geos/bin",
  "/opt/homebrew/opt/proj/bin",
  Sys.getenv("PATH"), sep = ":"))

cran_pkgs <- c(
  "yaml", "readxl", "dplyr", "tidyr", "purrr", "stringr",
  "ggplot2", "scales", "patchwork", "ggridges",
  "systemfonts", "ragg",
  "IntCal", "rcarbon", "posterior", "loo", "bayesplot",
  "cmdstanr"
)

# data.table built from source against Homebrew libomp fails to load on this
# machine. The CRAN macOS binary works. cmdstanr depends on it.
if (!requireNamespace("data.table", quietly = TRUE)) {
  try(install.packages("data.table", type = "binary"))
}

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
