# Bayesian Assessment of Contemporaneity in Radiocarbon Dates

A comprehensive R and Stan pipeline for analyzing radiocarbon dates from multiple archaeological deposits to assess temporal contemporaneity using Bayesian hierarchical models.

## Overview

This project provides a complete, documented framework for:

1. **Radiocarbon calibration** using IntCal20 via `rcarbon`
2. **Bayesian hierarchical modeling** in Stan to assess deposit contemporaneity
3. **Model comparison** using LOO cross-validation
4. **Simulation-based validation** with known ground truth
5. **Comprehensive visualization** of results and diagnostics

## Project Structure

```
project/
├── R/                          # R functions
│   ├── setup.R                 # Package management and setup
│   ├── data_validation.R       # Data validation functions
│   ├── calibration.R           # Radiocarbon calibration wrappers
│   ├── simulation.R            # Simulation framework
│   ├── stan_utilities.R        # Stan interface functions
│   ├── contemporaneity_metrics.R  # Metric calculation
│   ├── plotting.R              # Visualization functions
│   └── validation_framework.R  # Model validation system
├── stan/                       # Stan models
│   ├── contemporaneous_model.stan
│   ├── sequential_model.stan
│   └── partial_overlap_model.stan
├── tests/                      # Unit tests
│   └── testthat/
│       ├── test_data_validation.R
│       ├── test_calibration.R
│       └── test_simulation.R
├── data/                       # Data directory
├── docs/                       # Documentation
├── analysis.qmd                # Main analysis document (Quarto)
└── README.md                   # This file
```

## Installation

### Prerequisites

- R (≥ 4.0.0)
- RStudio (optional, recommended)
- C++ compiler (for Stan)

### Quick Start

```r
# Clone or download this repository
# Set working directory to project root

# Run setup script
source("R/setup.R")
complete_setup()
```

This will:
- Install all required R packages
- Install CmdStan toolchain
- Verify Stan installation
- Print package versions

### Required Packages

- `rcarbon` - Radiocarbon calibration
- `cmdstanr` - Stan interface (recommended over `rstan`)
- `posterior` - Posterior distribution analysis
- `bayesplot` - Diagnostic plots
- `loo` - Model comparison
- `dplyr`, `tidyr` - Data manipulation
- `ggplot2` - Visualization
- `testthat` - Unit testing
- `quarto` - Document rendering

## Usage

### Quick Example

```r
# Load functions
source("R/setup.R")
source("R/data_validation.R")
source("R/calibration.R")
source("R/simulation.R")
source("R/stan_utilities.R")
source("R/contemporaneity_metrics.R")
source("R/plotting.R")

# 1. Prepare data
c14_data <- data.frame(
  lab_code = c("OxA-1", "OxA-2", "Beta-1", "Beta-2"),
  age = c(2450, 2480, 2520, 2550),
  error = c(30, 35, 40, 30),
  deposit = c("A", "A", "B", "B")
)

# 2. Validate
validate_c14_data(c14_data)

# 3. Calibrate
cal_result <- calibrate_dates(
  ages = c14_data$age,
  errors = c14_data$error,
  lab_codes = c14_data$lab_code
)

# 4. Prepare for Stan
stan_data <- prepare_stan_data(c14_data, cal_result)

# 5. Compile and fit model
model <- compile_stan_model("stan/partial_overlap_model.stan")
fit <- fit_stan_model(model, stan_data)

# 6. Analyze results
metrics <- compute_all_metrics(fit, n_deposits = 2)
print_metrics_summary(metrics)

# 7. Visualize
plot_temporal_overlap(fit, n_deposits = 2)
plot_overlap_matrix(metrics$overlap_matrix)
```

### Running Full Analysis

To run the complete analysis workflow:

```bash
quarto render analysis.qmd
```

This generates an HTML report with all analyses, figures, and tables.

## Models

### 1. Contemporaneous Model

**Assumption**: All deposits share a single occupation window.

**Use case**: Testing hypothesis that all deposits are fully contemporaneous.

**Parameters**:
- `theta_start`: Shared start boundary (cal BP)
- `theta_end`: Shared end boundary (cal BP)
- `calendar_dates[n]`: True calendar date for each sample

### 2. Sequential Model

**Assumption**: Deposits are temporally ordered with no overlap.

**Use case**: Testing hypothesis of stratigraphic succession.

**Parameters**:
- `theta_start[k]`: Start boundary for deposit k
- `theta_end[k]`: End boundary for deposit k
- Constraint: `theta_end[k] < theta_start[k+1]`

### 3. Partial Overlap Model (Most Flexible)

**Assumption**: Each deposit has independent boundaries (allows any overlap pattern).

**Use case**: Default model for contemporaneity assessment.

**Parameters**:
- `theta_start[k]`: Start boundary for deposit k
- `theta_end[k]`: End boundary for deposit k
- No ordering constraints

## Model Comparison

Models are compared using PSIS-LOO cross-validation:

```r
comparison <- compare_models_loo(
  contemporaneous = fit_contemp,
  sequential = fit_seq,
  partial = fit_partial
)
```

**Interpretation**:
- **LOO weight > 0.9**: Strong evidence for model
- **LOO weight > 0.7**: Moderate evidence
- **LOO weight < 0.7**: Weak evidence; consider model averaging

## Validation

The framework includes comprehensive validation via simulation:

```r
# Run validation study
validation_results <- run_validation_study(
  n_trials = 100,
  scenarios = c("contemporaneous", "sequential", "partial"),
  sample_sizes = c(5, 10, 20),
  parallel = TRUE
)

# Summarize results
summary <- summarize_validation_results(validation_results)
print_validation_summary(summary)

# Visualize
plots <- plot_validation_study(validation_results)
print(plots$accuracy)
print(plots$coverage)
```

## Metrics

The following contemporaneity metrics are computed:

### Pairwise Overlap
- **Overlap probability**: P(deposits i and j overlap)
- **Overlap duration**: Expected length of overlap (years)

### Occupation Spans
- **Duration**: Length of occupation (with credible intervals)
- **Boundaries**: Start and end dates (with uncertainty)

### Multi-Deposit
- **Global overlap**: P(all deposits share common window)
- **Total time span**: Range encompassing all deposits

### Temporal Ordering
- **Sequence probability**: P(specific temporal order)

## Interpretation Guidelines

### Sample Size Recommendations

Based on validation studies:

| Scenario | Minimum Samples | Recommended |
|----------|----------------|-------------|
| Well-separated deposits | 5 per deposit | 10+ |
| Partially overlapping | 10 per deposit | 15+ |
| Fully contemporaneous | 8 per deposit | 12+ |

### Statistical vs. Cultural Contemporaneity

**Important**: Statistical overlap does not necessarily imply cultural contemporaneity:

- Statistical contemporaneity = temporal overlap within measurement uncertainty
- Cultural contemporaneity = actual simultaneous occupation or use
- Consider taphonomic processes, sample context, and archaeological interpretation

## Testing

Run unit tests:

```r
library(testthat)
test_dir("tests/testthat")
```

Tests cover:
- Data validation (edge cases, error checking)
- Simulation (parameter correctness, reproducibility)
- Calibration (numerical accuracy, structure)
- Stan utilities (integration)

## Reproducibility

For full reproducibility:

1. **Set random seeds** in all analysis scripts
2. **Record package versions**: `sessionInfo()`
3. **Document Stan version**: `cmdstanr::cmdstan_version()`
4. **Version control** all code and data
5. **Use relative file paths**

Example:

```r
set.seed(42)
sessionInfo()
cmdstanr::cmdstan_version()
```

## Calibration Curve

By default, IntCal20 is used [@reimer2020intcal20]. To use a different curve:

```r
cal_result <- calibrate_dates(
  ages = ages,
  errors = errors,
  calCurve = "intcal13"  # or "shcal20" for Southern Hemisphere
)
```

**Important**: IntCal20 is appropriate for Northern Hemisphere terrestrial samples. For marine samples or Southern Hemisphere, use appropriate curves and reservoir corrections.

## Troubleshooting

### Stan Compilation Errors

If Stan models fail to compile:

```r
# Reinstall CmdStan
cmdstanr::install_cmdstan(overwrite = TRUE)

# Check C++ compiler
cmdstanr::check_cmdstan_toolchain()
```

### Convergence Issues

If models show poor convergence:

1. Increase `adapt_delta`: `adapt_delta = 0.99`
2. Increase iterations: `iter_warmup = 2000, iter_sampling = 2000`
3. Increase `max_treedepth`: `max_treedepth = 15`
4. Check for label switching or multimodality
5. Provide better initial values

### Memory Issues

For large datasets:

- Use `method = "mixture"` for calibration (faster)
- Reduce number of posterior samples
- Use parallelization carefully

## Citation

If you use this framework, please cite:

```
[Project Citation]

And the key software packages:
- Stan: Carpenter et al. (2017)
- rcarbon: Crema & Bevan (2021)
- IntCal20: Reimer et al. (2020)
```

## References

1. Bronk Ramsey, C. (2009). Bayesian analysis of radiocarbon dates. *Radiocarbon*, 51(1), 337-360.

2. Carpenter, B., et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).

3. Crema, E. R., & Bevan, A. (2021). Inference from large sets of radiocarbon dates: software and methods. *Radiocarbon*, 63(1), 23-39.

4. Reimer, P. J., et al. (2020). The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0–55 cal kBP). *Radiocarbon*, 62(4), 725-757.

5. Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413-1432.

## License

[Specify license]

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Contact

[Contact information]

## Acknowledgments

This framework builds on methods developed by Christopher Bronk Ramsey (OxCal), the Stan Development Team, and the rcarbon package developers.
