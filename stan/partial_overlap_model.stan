/*
 * Model: Partial Overlap Model (Most Flexible)
 *
 * Description:
 *   Assumes each deposit has its own occupation window, but does NOT
 *   enforce any ordering constraints. This allows for any pattern of
 *   overlap between deposits (full, partial, or none).
 *
 * Parameters:
 *   - theta_start: Start dates for each deposit (cal BP) [n_deposits]
 *   - theta_end: End dates for each deposit (cal BP) [n_deposits]
 *   - calendar_dates: True calendar dates for each sample [N]
 *
 * Use:
 *   This is the most flexible model and should be used as the baseline
 *   for model comparison. It can represent any temporal relationship
 *   between deposits.
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
  int<lower=1> n_deposits;

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
  real<lower=0> prior_duration_mean;
  real<lower=0> prior_duration_sd;
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
  // Deposit-specific occupation boundaries (no ordering constraint)
  vector[n_deposits] theta_start;
  vector[n_deposits] theta_end;

  // True calendar dates for each sample
  vector[N] calendar_dates_raw;
}

transformed parameters {
  // Constrain calendar dates to deposit-specific windows
  vector[N] calendar_dates;
  vector<lower=0>[n_deposits] durations;

  // Calculate durations
  for (k in 1:n_deposits) {
    durations[k] = theta_start[k] - theta_end[k];
  }

  // Transform calendar dates to be within deposit windows
  for (n in 1:N) {
    int dep = deposit_id[n];
    // Transform raw parameter to be within [theta_end[dep], theta_start[dep]]
    calendar_dates[n] = theta_end[dep] +
                        durations[dep] * inv_logit(calendar_dates_raw[n]);
  }
}

model {
  // Priors on boundaries
  theta_start ~ normal(prior_start_mean, prior_start_sd);

  // Priors on durations (indirectly through theta_end)
  durations ~ normal(prior_duration_mean, prior_duration_sd);

  // Implicit prior: theta_start > theta_end (enforced through durations)
  for (k in 1:n_deposits) {
    target += -log(durations[k] + epsilon);  // Jacobian for positivity constraint
  }

  // Prior on calendar dates
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

  // Overlap probabilities between all pairs of deposits
  array[n_deposits, n_deposits] real overlap_prob;

  // Overlap durations
  array[n_deposits, n_deposits] real overlap_duration;

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

  // Compute pairwise overlaps
  for (i in 1:n_deposits) {
    for (j in 1:n_deposits) {
      if (i == j) {
        overlap_prob[i, j] = 1.0;
        overlap_duration[i, j] = durations[i];
      } else {
        // Check if deposits overlap
        // Deposits overlap if: theta_end[i] < theta_start[j] AND theta_end[j] < theta_start[i]
        // Equivalently: max(theta_end) < min(theta_start)
        real max_end = max(theta_end[i], theta_end[j]);
        real min_start = min(theta_start[i], theta_start[j]);

        if (max_end < min_start) {
          // There is overlap
          overlap_prob[i, j] = 1.0;
          overlap_duration[i, j] = min_start - max_end;
        } else {
          // No overlap
          overlap_prob[i, j] = 0.0;
          overlap_duration[i, j] = 0.0;
        }
      }
    }
  }
}
