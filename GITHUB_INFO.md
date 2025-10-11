# GitHub Repository Information

## Repository
**https://github.com/clipo/village**

## Successfully Pushed

✅ **Initial commit**: Complete radiocarbon dating contemporaneity analysis framework
- All R functions (9 files)
- All Stan models (3 files)
- All tests (3 files)
- All documentation (7+ files)
- All analysis scripts (4 files)

✅ **Second commit**: Sample radiocarbon dataset
- data/radiocarbon_dates.csv (157 dates from 47 sites)

## Repository Contents

```
village/
├── R/                          # 9 R function files (~92 KB)
├── stan/                       # 3 Stan models (~14 KB)
├── tests/testthat/             # 3 test files (~12 KB)
├── data/                       # Sample data (1 CSV file)
├── output/                     # Output directory (empty)
├── docs/                       # Additional docs directory
├── Documentation files (7)
├── Analysis scripts (4)
└── Configuration files
```

## Commits

1. **3e44eb6** - Initial commit with complete framework
2. **44c5547** - Add sample radiocarbon dataset

## Branch

- **main** (default branch)

## Total Files Pushed

29 files including:
- 9 R function modules
- 3 Stan models
- 3 test files
- 7 documentation files
- 4 analysis scripts
- 1 data file
- 2 configuration files

## How to Clone

```bash
git clone https://github.com/clipo/village.git
cd village
```

## Quick Start After Cloning

```r
# Verify installation
source("verify_installation.R")

# Preview your data
source("preview_data.R")

# Run analysis
source("analyze_your_data.R")
```

## Repository Structure Highlights

### Documentation
- **README.md** - Main documentation
- **QUICKSTART.md** - Getting started in 5 minutes
- **DATA_LOADING_GUIDE.md** - How to load your data
- **PROJECT_SUMMARY.md** - Technical overview
- **DELIVERABLES.md** - Complete inventory

### Key Scripts
- **preview_data.R** - Explore your data
- **analyze_your_data.R** - Complete analysis pipeline
- **example_workflow.R** - Examples with simulated data
- **verify_installation.R** - Check setup

### Core Functions
- **R/data_loading.R** - Load CSV/Excel files
- **R/calibration.R** - Radiocarbon calibration
- **R/stan_utilities.R** - Stan model interface
- **R/contemporaneity_metrics.R** - Calculate overlap
- **R/plotting.R** - Visualizations

### Models
- **stan/contemporaneous_model.stan** - Shared window
- **stan/sequential_model.stan** - Ordered deposits
- **stan/partial_overlap_model.stan** - Flexible model

## Features

✓ Load radiocarbon data from CSV/Excel
✓ Automatic column detection and mapping
✓ Calibrate with IntCal20 via rcarbon
✓ Three Bayesian hierarchical models in Stan
✓ Model comparison via LOO cross-validation
✓ Contemporaneity metrics and probabilities
✓ Publication-quality visualizations
✓ Complete validation framework
✓ Unit tests for quality assurance
✓ Comprehensive documentation

## Statistics

- **Total code**: ~177 KB
- **Lines**: ~5,180
- **Functions**: 53+ documented
- **Models**: 3 Stan programs
- **Tests**: 30+ unit tests
- **Documentation**: 850+ lines

## License

[To be specified]

## Citation

If you use this framework, please cite:
- This repository
- Stan (Carpenter et al. 2017)
- rcarbon (Crema & Bevan 2021)
- IntCal20 (Reimer et al. 2020)

## Contact

Repository: https://github.com/clipo/village
Issues: https://github.com/clipo/village/issues

## Last Updated

2025-10-11

---

**Status**: ✅ Complete and ready to use
**Branch**: main
**Commits**: 2
**Files**: 29
