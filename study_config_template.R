# Study Configuration Template
# Copy this file and customize it for your study
# This information will be incorporated into the generated paper

# =============================================================================
# STUDY METADATA
# =============================================================================

study_title <- "YOUR STUDY TITLE HERE"
study_subtitle <- "A Bayesian Approach to Radiocarbon Chronology"

author_name <- "Your Name"
author_affiliation <- "Your Institution/Department"
author_email <- "your.email@institution.edu"

# =============================================================================
# ARCHAEOLOGICAL CONTEXT
# =============================================================================

# Study region and period
study_region <- "Your study region (e.g., Thames Valley, southern England)"
study_period <- "Time period (e.g., Middle Bronze Age, ca. 3500-3000 BP)"
study_culture <- "Cultural context (e.g., Bronze Age settlement complexes)"

# Research questions
# Provide 2-4 specific research questions as a character vector
research_questions <- c(
  "Were the sites occupied simultaneously?",
  "Do sites cluster into distinct chronological phases?",
  "What is the duration of occupation at each site?",
  "How do occupation patterns relate to environmental/cultural changes?"
)

# =============================================================================
# SITE DESCRIPTIONS
# =============================================================================

# Provide a named list with descriptions for each site
# Keys should match site names in your data file
site_descriptions <- list(
  "Site_Name_1" = "Brief description of Site 1 (location, size, features, materials recovered, etc.)",
  "Site_Name_2" = "Brief description of Site 2",
  "Site_Name_3" = "Brief description of Site 3",
  "Site_Name_4" = "Brief description of Site 4",
  "Site_Name_5" = "Brief description of Site 5"
  # Add more as needed
)

# =============================================================================
# ARCHAEOLOGICAL CONTEXT (BACKGROUND)
# =============================================================================

# Provide 1-3 paragraphs describing the archaeological background
# This will be inserted into the Introduction section
archaeological_context <- "
[Paragraph 1: Regional/temporal context]
Describe the broader archaeological period and region. What is known about
this time period? What are the key cultural characteristics? What major
archaeological questions or debates exist?

[Paragraph 2: Site context]
Describe the specific sites in your study. How were they discovered? What
kinds of features and materials are present? What previous work has been done?

[Paragraph 3: Importance]
Explain why understanding contemporaneity is important for this case. What
would it mean if sites were contemporaneous vs. sequential? How does this
relate to broader research questions about the period/region?
"

# =============================================================================
# INTERPRETATION GUIDELINES
# =============================================================================

# What would contemporaneity mean for your research questions?
contemporaneity_implications <- "
If these sites were contemporaneous, it would suggest:
1. [Interpretation 1 - e.g., regional settlement network]
2. [Interpretation 2 - e.g., sufficient resources for multiple communities]
3. [Interpretation 3 - e.g., social complexity/hierarchy]
4. [Interpretation 4 - e.g., cultural interaction/exchange]

Temporal separation would instead indicate:
1. [Alternative interpretation 1 - e.g., sequential occupation phases]
2. [Alternative interpretation 2 - e.g., population movement]
3. [Alternative interpretation 3 - e.g., changing settlement strategies]
"

# =============================================================================
# ADDITIONAL DISCUSSION POINTS
# =============================================================================

# Optional: Additional contextual information for the Discussion section
# Provide as a named list
discussion_points <- list(
  material_culture = "Brief description of material culture patterns (ceramics, lithics, etc.)",
  subsistence = "Brief description of subsistence evidence (floral/faunal remains, etc.)",
  landscape = "Brief description of landscape/environmental context",
  regional_context = "How do these sites fit into broader regional patterns?",
  previous_chronology = "What was the previous understanding of site chronology?"
)

# =============================================================================
# DATA SPECIFICATIONS (OPTIONAL)
# =============================================================================

# If you want to override defaults for the analysis
analysis_options <- list(
  min_dates_per_site = 3,           # Minimum dates required per site
  mcmc_iterations_warmup = 1000,    # MCMC warmup iterations
  mcmc_iterations_sampling = 2000,  # MCMC sampling iterations
  mcmc_chains = 4,                  # Number of MCMC chains
  calibration_curve = "intcal20",   # Calibration curve (intcal20, marine20, shcal20)
  figure_dpi = 300,                 # DPI for saved figures
  output_directory = "output/my_study"  # Where to save results
)

# =============================================================================
# NOTES
# =============================================================================

# After customizing this file:
# 1. Save with a descriptive name (e.g., my_bronze_age_study_config.R)
# 2. Source it in your analysis script: source("my_bronze_age_study_config.R")
# 3. Use it to generate a custom paper: generate_custom_paper(config_file = "my_bronze_age_study_config.R", ...)
