# Quick Start Guide

Get up and running with radiocarbon dating contemporaneity analysis in 5 minutes.

## Step 1: Verify Installation (2 minutes)

```r
# Run verification script
source("verify_installation.R")
```

If checks fail, install missing components:

```r
# Install R packages
source("R/setup.R")
complete_setup()

# Install CmdStan (if needed)
library(cmdstanr)
install_cmdstan()
```

## Step 2: Run Example Analysis (3 minutes)

### Option A: Interactive Analysis

Open R or RStudio and run:

```r
# Load functions
source("R/setup.R")
source("R/data_validation.R")
source("R/calibration.R")
source("R/simulation.R")
source("R/stan_utilities.R")
source("R/contemporaneity_metrics.R")
source("R/plotting.R")

# Generate example data
sim <- simulate_partial_overlap(
  n_deposits = 2,
  overlap_percentage = 60,
  deposit_duration = 250,
  earliest_start = 2800,
  n_dates_per_deposit = 10,
  measurement_error = 30,
  seed = 42
)

# Extract radiocarbon data
c14_data <- sim$dates[, c("lab_code", "age", "error", "deposit")]

# Calibrate
cal_result <- calibrate_dates(c14_data$age, c14_data$error, c14_data$lab_code)

# Prepare for Stan
stan_data <- prepare_stan_data(c14_data, cal_result)

# Compile and fit model
model <- compile_stan_model("stan/partial_overlap_model.stan")
fit <- fit_stan_model(model, stan_data, chains = 2, iter_sampling = 500)

# Analyze
metrics <- compute_all_metrics(fit, n_deposits = 2)
print_metrics_summary(metrics)

# Visualize
plot_temporal_overlap(fit, n_deposits = 2)
```

### Option B: Complete Workflow Script

```r
# Run the full example workflow
source("example_workflow.R")
```

This will:
- Generate simulated data
- Fit all three models
- Compare models
- Create all plots
- Save outputs to `output/` directory

### Option C: Generate Full Report

```bash
quarto render analysis.qmd
```

This creates a complete HTML report with all analyses.

## Step 3: Analyze Your Own Data

Create your data frame:

```r
my_data <- data.frame(
  lab_code = c("OxA-12345", "OxA-12346", "Beta-54321"),
  age = c(2450, 2480, 2520),
  error = c(30, 35, 40),
  deposit = c("Deposit_A", "Deposit_A", "Deposit_B")
)

# Validate
validate_c14_data(my_data)

# Calibrate
cal_result <- calibrate_dates(
  ages = my_data$age,
  errors = my_data$error,
  lab_codes = my_data$lab_code
)

# Continue with workflow...
```

## Common Tasks

### View Calibrated Dates

```r
plot_calibrated_dates(cal_result, my_data)
```

### Compare Multiple Models

```r
# Fit all models
fit_contemp <- fit_stan_model(model_contemp, stan_data)
fit_seq <- fit_stan_model(model_seq, stan_data)
fit_partial <- fit_stan_model(model_partial, stan_data)

# Compare
comparison <- compare_models_loo(
  contemporaneous = fit_contemp,
  sequential = fit_seq,
  partial = fit_partial
)

# Best model
best <- names(comparison$weights)[which.max(comparison$weights)]
print(paste("Best model:", best))
```

### Calculate Overlap Probability

```r
# For two deposits
prob <- calculate_overlap_probability(fit, deposit_i = 1, deposit_j = 2)
print(paste("Overlap probability:", round(prob, 3)))
```

### Create Overlap Matrix

```r
overlap_matrix <- create_overlap_matrix(fit, n_deposits = 3)
plot_overlap_matrix(overlap_matrix)
```

### Check Convergence

```r
check_stan_diagnostics(fit)
```

### Posterior Predictive Check

```r
plot_posterior_predictive(fit, my_data)
```

## Troubleshooting

### Stan Won't Compile

```r
# Check toolchain
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Reinstall CmdStan
cmdstanr::install_cmdstan(overwrite = TRUE)
```

### Models Don't Converge

```r
# Increase adapt_delta
fit <- fit_stan_model(
  model,
  stan_data,
  adapt_delta = 0.99,  # Default is 0.95
  iter_warmup = 2000   # Default is 1000
)
```

### Out of Memory

```r
# Use fewer posterior samples
fit <- fit_stan_model(
  model,
  stan_data,
  chains = 2,          # Instead of 4
  iter_sampling = 500  # Instead of 1000
)
```

### Package Installation Issues

```r
# Specify CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install one at a time
install.packages("rcarbon")
install.packages("cmdstanr",
                 repos = c("https://mc-stan.org/r-packages/",
                          getOption("repos")))
```

## Next Steps

1. **Read Documentation**: `README.md` has complete details
2. **Review Examples**: `example_workflow.R` shows full pipeline
3. **Read Methods**: `analysis.qmd` explains the approach
4. **Run Validation**: Test model performance on your data scale

## Data Requirements

**Minimum**:
- 5 dates per deposit
- 2 deposits
- Standard measurement errors (15-50 years typical)

**Recommended**:
- 10-15 dates per deposit
- Known material type
- Stratigraphic context
- Understand potential taphonomic issues

## Interpretation Checklist

Before concluding deposits are contemporaneous:

- [ ] Models converged (Rhat < 1.01, ESS > 400)
- [ ] Posterior predictive check looks good
- [ ] LOO comparison favors a specific model
- [ ] Overlap probability is high (>0.7 for "likely")
- [ ] Credible intervals are reasonable widths
- [ ] Archaeological context supports conclusion
- [ ] Taphonomic issues considered
- [ ] Reservoir effects ruled out (if applicable)

## Getting Help

1. Check `README.md` troubleshooting section
2. Review `PROJECT_SUMMARY.md` for overview
3. Examine `example_workflow.R` for working code
4. Run `verify_installation.R` to diagnose setup
5. Check function documentation: `?function_name`
6. Review Stan model diagnostics carefully

## Key Functions Reference

| Function | Purpose |
|----------|---------|
| `validate_c14_data()` | Check data validity |
| `calibrate_dates()` | Calibrate radiocarbon dates |
| `compile_stan_model()` | Compile Stan program |
| `fit_stan_model()` | Fit Bayesian model |
| `compare_models_loo()` | Model comparison |
| `compute_all_metrics()` | Calculate contemporaneity metrics |
| `plot_temporal_overlap()` | Visualize occupation ranges |
| `create_overlap_matrix()` | Pairwise overlap probabilities |

## File Locations

- **R functions**: `R/*.R`
- **Stan models**: `stan/*.stan`
- **Tests**: `tests/testthat/*.R`
- **Examples**: `example_workflow.R`
- **Documentation**: `README.md`, `analysis.qmd`
- **Output**: `output/` (created automatically)

## Session Info

Always record at the end of your analysis:

```r
sessionInfo()
cmdstanr::cmdstan_version()
```

---

**Ready to start?** Run `source("verify_installation.R")` and then `source("example_workflow.R")`!
