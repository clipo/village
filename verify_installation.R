#' Verify Installation and Setup
#'
#' This script verifies that all components are properly installed
#' and ready for analysis.
#'
#' @author Project Team
#' @date 2025-10-11

cat("==========================================================\n")
cat("  Radiocarbon Dating Analysis Framework\n")
cat("  Installation Verification\n")
cat("==========================================================\n\n")

# Track success
all_checks_passed <- TRUE

# ==============================================================================
# 1. Check R Version
# ==============================================================================

cat("1. Checking R version...\n")
r_version <- R.version$version.string
cat("   ", r_version, "\n")

if (getRversion() < "4.0.0") {
  cat("   WARNING: R >= 4.0.0 recommended\n")
  all_checks_passed <- FALSE
} else {
  cat("   ✓ R version OK\n")
}

# ==============================================================================
# 2. Check Required Packages
# ==============================================================================

cat("\n2. Checking required R packages...\n")

required_packages <- c(
  "rcarbon", "cmdstanr", "posterior", "bayesplot", "loo",
  "dplyr", "tidyr", "ggplot2", "testthat", "quarto"
)

missing_packages <- character(0)

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- as.character(packageVersion(pkg))
    cat(sprintf("   ✓ %-15s : %s\n", pkg, ver))
  } else {
    cat(sprintf("   ✗ %-15s : NOT INSTALLED\n", pkg))
    missing_packages <- c(missing_packages, pkg)
    all_checks_passed <- FALSE
  }
}

if (length(missing_packages) > 0) {
  cat("\n   To install missing packages, run:\n")
  cat("   install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n", sep = "")
}

# ==============================================================================
# 3. Check Stan Installation
# ==============================================================================

cat("\n3. Checking Stan installation...\n")

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  library(cmdstanr)

  # Check if CmdStan is installed
  cmdstan_installed <- tryCatch({
    ver <- cmdstan_version()
    cat("   ✓ CmdStan version: ", ver, "\n")
    TRUE
  }, error = function(e) {
    cat("   ✗ CmdStan not found\n")
    cat("     To install: cmdstanr::install_cmdstan()\n")
    FALSE
  })

  if (!cmdstan_installed) {
    all_checks_passed <- FALSE
  } else {
    # Check C++ toolchain
    toolchain_ok <- tryCatch({
      check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
      cat("   ✓ C++ toolchain OK\n")
      TRUE
    }, error = function(e) {
      cat("   ✗ C++ toolchain issue\n")
      cat("     Run: cmdstanr::check_cmdstan_toolchain()\n")
      FALSE
    })

    if (!toolchain_ok) {
      all_checks_passed <- FALSE
    }
  }
} else {
  cat("   ✗ cmdstanr not installed\n")
  all_checks_passed <- FALSE
}

# ==============================================================================
# 4. Check Project Structure
# ==============================================================================

cat("\n4. Checking project structure...\n")

required_dirs <- c("R", "stan", "tests", "data", "docs", "output")
required_files <- c(
  "R/setup.R",
  "R/data_validation.R",
  "R/calibration.R",
  "R/simulation.R",
  "R/stan_utilities.R",
  "R/contemporaneity_metrics.R",
  "R/plotting.R",
  "R/validation_framework.R",
  "stan/contemporaneous_model.stan",
  "stan/sequential_model.stan",
  "stan/partial_overlap_model.stan",
  "analysis.qmd",
  "example_workflow.R",
  "README.md"
)

# Check directories
for (dir in required_dirs) {
  if (dir.exists(dir)) {
    cat(sprintf("   ✓ %s/\n", dir))
  } else {
    cat(sprintf("   ✗ %s/ NOT FOUND\n", dir))
    all_checks_passed <- FALSE
  }
}

# Check key files
for (file in required_files) {
  if (file.exists(file)) {
    cat(sprintf("   ✓ %s\n", file))
  } else {
    cat(sprintf("   ✗ %s NOT FOUND\n", file))
    all_checks_passed <- FALSE
  }
}

# ==============================================================================
# 5. Test Basic Functionality
# ==============================================================================

cat("\n5. Testing basic functionality...\n")

# Test data validation
test_validation <- tryCatch({
  source("R/data_validation.R", local = TRUE)
  test_data <- data.frame(
    lab_code = c("TEST-1", "TEST-2"),
    age = c(2450, 2480),
    error = c(30, 35),
    deposit = c("A", "A")
  )
  validate_c14_data(test_data)
  cat("   ✓ Data validation works\n")
  TRUE
}, error = function(e) {
  cat("   ✗ Data validation error:", e$message, "\n")
  FALSE
})

if (!test_validation) {
  all_checks_passed <- FALSE
}

# Test calibration (requires rcarbon)
test_calibration <- tryCatch({
  source("R/calibration.R", local = TRUE)
  cal_curve <- get_calibration_curve("intcal20")
  cat("   ✓ Calibration functions work\n")
  cat("     IntCal20 curve loaded:", nrow(cal_curve), "points\n")
  TRUE
}, error = function(e) {
  cat("   ✗ Calibration error:", e$message, "\n")
  FALSE
})

if (!test_calibration) {
  all_checks_passed <- FALSE
}

# Test Stan model compilation (if Stan available)
if (requireNamespace("cmdstanr", quietly = TRUE) && cmdstan_installed) {
  test_stan <- tryCatch({
    source("R/stan_utilities.R", local = TRUE)
    model <- cmdstan_model("stan/contemporaneous_model.stan", compile = FALSE)
    cat("   ✓ Stan models can be loaded\n")
    TRUE
  }, error = function(e) {
    cat("   ✗ Stan model error:", e$message, "\n")
    FALSE
  })

  if (!test_stan) {
    all_checks_passed <- FALSE
  }
}

# ==============================================================================
# 6. Summary
# ==============================================================================

cat("\n==========================================================\n")
if (all_checks_passed) {
  cat("  ✓ ALL CHECKS PASSED\n")
  cat("  Installation is complete and ready to use!\n")
  cat("\n")
  cat("  Next steps:\n")
  cat("  1. Review README.md for usage instructions\n")
  cat("  2. Run example_workflow.R to test the pipeline\n")
  cat("  3. Or render analysis.qmd for full demonstration\n")
} else {
  cat("  ✗ SOME CHECKS FAILED\n")
  cat("  Please address the issues above before proceeding.\n")
  cat("\n")
  cat("  Common solutions:\n")
  cat("  - Install missing packages: source('R/setup.R'); complete_setup()\n")
  cat("  - Install CmdStan: cmdstanr::install_cmdstan()\n")
  cat("  - Check R version (>= 4.0.0 recommended)\n")
}
cat("==========================================================\n")

# Return status invisibly
invisible(all_checks_passed)
