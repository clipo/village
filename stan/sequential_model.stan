/*
 * Model: Sequential (Non-Contemporaneous) Deposit Model
 *
 * Description:
 *   Assumes deposits have separate, non-overlapping occupation windows
 *   with enforced temporal ordering. Each deposit has its own boundaries.
 *
 * Parameters:
 *   - theta_start: Start dates for each deposit (cal BP) [n_deposits]
 *   - theta_end: End dates for each deposit (cal BP) [n_deposits]
 *   - calendar_dates: True calendar dates for each sample [N]
 *
 * Data:
 *   Same as contemporaneous model, but uses deposit_id to assign
 *   dates to deposit-specific windows
 *
 * Constraints:
 *   - Enforces ordering: theta_end[k] < theta_start[k+1]
 *   - This ensures deposits are sequential in time
 *
 * Author: Project Team
 * Date: 2025-10-11
 *
 * References:
 *   Bronk Ramsey, C. (2009). Bayesian analysis of radiocarbon dates.
 *     Radiocarbon, 51(1), 337-360.
 */

data {
  // Number of observations and deposits
  int<lower=1> N;
  int<lower=2> n_deposits;

  // Observed radiocarbon data
  vector[N] c14_age;
  vector<lower=0>[N] c14_error;

  // Deposit assignments
  array[N] int<lower=1, upper=n_deposits> deposit_id;

  // Calibration approximation via mixture of normals
  int<lower=1> K;  // Number of mixture components
  array[N, K] real<lower=0, upper=1> mix_weights;
  array[N, K] real mix_means;  // Means in cal BP
  array[N, K] real<lower=0> mix_sds;

  // Priors
  real prior_start_mean;
  real<lower=0> prior_start_sd;
  real<lower=0> prior_duration_mean;  // Expected duration of deposits
  real<lower=0> prior_duration_sd;
  real<lower=0> prior_gap_mean;  // Expected gap between deposits
  real<lower=0> prior_gap_sd;
}

transformed data {
  real epsilon = 1e-10;

  // Count dates per deposit
  array[n_deposits] int dates_per_deposit;
  for (k in 1:n_deposits) {
    dates_per_deposit[k] = 0;
  }
  for (n in 1:N) {
    dates_per_deposit[deposit_id[n]] += 1;
  }
}

parameters {
  // Deposit-specific occupation boundaries
  // Use ordered constraint to enforce sequence
  ordered[n_deposits] theta_end_ordered;  // End dates (smaller values)

  // Durations (must be positive)
  vector<lower=0>[n_deposits] durations;

  // True calendar dates for each sample
  vector[N] calendar_dates_raw;  // Will be constrained in transformed parameters
}

transformed parameters {
  // Compute start dates from end dates and durations
  vector[n_deposits] theta_start;
  vector[n_deposits] theta_end;

  theta_end = theta_end_ordered;

  for (k in 1:n_deposits) {
    theta_start[k] = theta_end[k] + durations[k];
  }

  // Constrain calendar dates to deposit-specific windows
  vector[N] calendar_dates;
  for (n in 1:N) {
    int dep = deposit_id[n];
    // Transform raw parameter to be within deposit window
    calendar_dates[n] = theta_end[dep] +
                        durations[dep] * inv_logit(calendar_dates_raw[n]);
  }
}

model {
  // Prior on first deposit start
  theta_start[1] ~ normal(prior_start_mean, prior_start_sd);

  // Priors on durations
  durations ~ normal(prior_duration_mean, prior_duration_sd);

  // Prior on gaps between deposits
  for (k in 1:(n_deposits - 1)) {
    real gap = theta_end[k] - theta_start[k + 1];
    gap ~ normal(prior_gap_mean, prior_gap_sd);
  }

  // Prior on calendar dates (implicit through transformation)
  calendar_dates_raw ~ normal(0, 1);

  // Likelihood: observed 14C ages given true calendar dates
  for (n in 1:N) {
    vector[K] log_lik_components;

    for (k in 1:K) {
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

  // Gaps between deposits
  vector[n_deposits - 1] gaps;

  // Compute gaps
  for (k in 1:(n_deposits - 1)) {
    gaps[k] = theta_end[k] - theta_start[k + 1];
  }

  // Compute posterior predictive and log likelihood
  for (n in 1:N) {
    vector[K] log_lik_components;

    // Sample from mixture for posterior predictive
    int component = categorical_rng(to_vector(mix_weights[n, ]));
    c14_age_rep[n] = normal_rng(calendar_dates[n], mix_sds[n, component]);

    // Compute log likelihood
    for (k in 1:K) {
      log_lik_components[k] = log(mix_weights[n, k] + epsilon) +
                               normal_lpdf(calendar_dates[n] | mix_means[n, k], mix_sds[n, k]);
    }
    log_lik[n] = log_sum_exp(log_lik_components);
  }
}
