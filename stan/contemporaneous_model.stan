/*
 * Model: Fully Contemporaneous Occupation Model
 *
 * Description:
 *   Assumes all deposits share a single occupation window. All radiocarbon
 *   dates from all deposits are drawn from the same temporal interval.
 *
 * Parameters:
 *   - theta_start: Start of shared occupation (cal BP) - earlier/higher value
 *   - theta_end: End of shared occupation (cal BP) - later/lower value
 *   - calendar_dates: True calendar dates for each sample (latent variables)
 *
 * Data:
 *   - N: Number of radiocarbon dates
 *   - c14_age: Observed radiocarbon ages (BP)
 *   - c14_error: Measurement errors (1-sigma)
 *   - deposit_id: Deposit assignment for each date (not used in this model)
 *   - K: Number of mixture components (for calibration approximation)
 *   - mix_weights: Mixture weights [N x K]
 *   - mix_means: Mixture means in cal BP [N x K]
 *   - mix_sds: Mixture standard deviations [N x K]
 *
 * Approach:
 *   Uses mixture of normals approximation to calibrated distributions.
 *   This is faster than including full calibration curve in Stan.
 *
 * Author: Project Team
 * Date: 2025-10-11
 *
 * References:
 *   Bronk Ramsey, C. (2009). Bayesian analysis of radiocarbon dates.
 *     Radiocarbon, 51(1), 337-360.
 */

data {
  // Number of observations
  int<lower=1> N;

  // Observed radiocarbon data
  vector[N] c14_age;
  vector<lower=0>[N] c14_error;

  // Deposit assignments (informational only for this model)
  array[N] int<lower=1> deposit_id;
  int<lower=1> n_deposits;

  // Calibration approximation via mixture of normals
  int<lower=1> K;  // Number of mixture components
  array[N, K] real<lower=0, upper=1> mix_weights;  // Must sum to 1 per row
  array[N, K] real mix_means;  // Means in cal BP
  array[N, K] real<lower=0> mix_sds;  // Standard deviations

  // Priors (optional, can use defaults)
  real prior_start_mean;  // Expected start of occupation
  real<lower=0> prior_start_sd;  // Uncertainty in start
  real prior_end_mean;  // Expected end of occupation
  real<lower=0> prior_end_sd;  // Uncertainty in end
}

transformed data {
  // For numerical stability
  real epsilon = 1e-10;
}

parameters {
  // Shared occupation boundaries
  real theta_start;  // Start of occupation (cal BP, larger value)
  real theta_end;    // End of occupation (cal BP, smaller value)

  // True calendar dates for each sample (latent variables)
  vector<lower=theta_end, upper=theta_start>[N] calendar_dates;
}

transformed parameters {
  // Duration of occupation
  real<lower=0> occupation_duration = theta_start - theta_end;
}

model {
  // Priors on occupation boundaries
  theta_start ~ normal(prior_start_mean, prior_start_sd);
  theta_end ~ normal(prior_end_mean, prior_end_sd);

  // Constraint: start must be before end (in cal BP convention)
  target += -log(theta_start - theta_end + epsilon);  // Jacobian adjustment

  // Prior on calendar dates: uniform within occupation window
  // (implicit through constraints)

  // Likelihood: observed 14C ages given true calendar dates
  for (n in 1:N) {
    // Use mixture approximation to calibrated distribution
    vector[K] log_lik_components;

    for (k in 1:K) {
      // Likelihood that true calendar date produces observed 14C age
      // via the calibration curve (approximated by mixture component)
      log_lik_components[k] = log(mix_weights[n, k] + epsilon) +
                               normal_lpdf(calendar_dates[n] | mix_means[n, k], mix_sds[n, k]);
    }

    target += log_sum_exp(log_lik_components);
  }
}

generated quantities {
  // Posterior predictive checks
  array[N] real c14_age_rep;

  // Log likelihood for LOO
  vector[N] log_lik;

  // Compute posterior predictive and log likelihood
  for (n in 1:N) {
    vector[K] log_lik_components;

    // Sample from mixture for posterior predictive
    int component = categorical_rng(to_vector(mix_weights[n, ]));
    c14_age_rep[n] = normal_rng(calendar_dates[n], mix_sds[n, component]);

    // Compute log likelihood for this observation
    for (k in 1:K) {
      log_lik_components[k] = log(mix_weights[n, k] + epsilon) +
                               normal_lpdf(calendar_dates[n] | mix_means[n, k], mix_sds[n, k]);
    }
    log_lik[n] = log_sum_exp(log_lik_components);
  }
}
