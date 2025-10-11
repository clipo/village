#' Setup and Package Management for Radiocarbon Dating Analysis
#'
#' This script manages package dependencies and verifies Stan installation
#' for the radiocarbon dating contemporaneity analysis pipeline.
#'
#' @author Project Team
#' @date 2025-10-11

# Required packages
required_packages <- c(
  "rcarbon",       # Radiocarbon calibration
  "cmdstanr",      # Stan interface (preferred for latest features)
  "posterior",     # Working with posterior distributions
  "bayesplot",     # Posterior visualization and diagnostics
  "loo",           # Model comparison via PSIS-LOO
  "dplyr",         # Data manipulation
  "tidyr",         # Data tidying
  "ggplot2",       # Visualization
  "testthat",      # Unit testing
  "quarto"         # Document generation
)

#' Install missing packages
#'
#' @param packages Character vector of package names
#' @return NULL (installs packages as side effect)
install_missing_packages <- function(packages) {
  installed <- installed.packages()[, "Package"]
  missing <- packages[!packages %in% installed]

  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  } else {
    message("All required packages are already installed.")
  }
}

#' Verify Stan installation
#'
#' Checks that cmdstanr is installed and CmdStan toolchain is available
#'
#' @return List with Stan version info and status
verify_stan_installation <- function() {
  if (!require("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is not installed. Install with: install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))")
  }

  # Check if CmdStan is installed
  cmdstan_installed <- tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!cmdstan_installed) {
    message("CmdStan toolchain not found. Installing...")
    cmdstanr::install_cmdstan()
  }

  # Get version information
  stan_info <- list(
    cmdstan_version = cmdstanr::cmdstan_version(),
    cmdstan_path = cmdstanr::cmdstan_path(),
    cmdstanr_version = as.character(packageVersion("cmdstanr"))
  )

  message("Stan installation verified:")
  message("  CmdStan version: ", stan_info$cmdstan_version)
  message("  CmdStan path: ", stan_info$cmdstan_path)
  message("  cmdstanr version: ", stan_info$cmdstanr_version)

  return(invisible(stan_info))
}

#' Get complete session info including all package versions
#'
#' @return Session info object with package versions
get_package_versions <- function() {
  sess_info <- sessionInfo()

  message("\n=== R Session Information ===")
  message("R version: ", R.version.string)
  message("Platform: ", sess_info$platform)
  message("Running under: ", sess_info$running)

  message("\n=== Package Versions ===")

  # Get versions of required packages
  for (pkg in required_packages) {
    if (pkg %in% installed.packages()[, "Package"]) {
      ver <- as.character(packageVersion(pkg))
      message(sprintf("  %-15s: %s", pkg, ver))
    } else {
      message(sprintf("  %-15s: NOT INSTALLED", pkg))
    }
  }

  return(invisible(sess_info))
}

#' Complete setup verification
#'
#' Runs all setup checks and reports status
#'
#' @return List with setup status
complete_setup <- function() {
  message("=== Radiocarbon Dating Analysis Setup ===\n")

  # Install missing packages
  install_missing_packages(required_packages)

  # Verify Stan
  stan_info <- verify_stan_installation()

  # Get package versions
  session_info <- get_package_versions()

  message("\n=== Setup Complete ===")
  message("All dependencies are installed and ready.")
  message("Run get_package_versions() anytime to check versions.")

  return(invisible(list(
    stan_info = stan_info,
    session_info = session_info
  )))
}

# Justification for cmdstanr over rstan:
# cmdstanr is chosen because it:
# 1. Provides access to the latest Stan features and bug fixes
# 2. Has better performance and more efficient memory usage
# 3. Simplifies installation and maintenance (no compilation during R package install)
# 4. Offers better support for modern Stan workflows (parallel chains, within-chain parallelization)
# 5. Provides cleaner interface for diagnostics and posterior analysis
# 6. Is the recommended interface for new Stan projects (as of 2024)

# If you need to use rstan instead (e.g., for compatibility reasons),
# replace cmdstanr with rstan in the required_packages and modify
# verify_stan_installation() accordingly.
