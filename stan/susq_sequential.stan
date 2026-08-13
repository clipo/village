/**
 * Strictly serial occupation. Every site has one phase and no two phases
 * overlap. Boundaries form a single ordered vector, so phase k ends before
 * phase k+1 begins.
 *
 * The order of sites is fixed by the data, ascending in median calibrated
 * date. That is the ordering most favourable to the serial hypothesis, chosen
 * using the data, so a win for this model is suggestive rather than decisive
 * while a loss is correspondingly strong. See spec section 5.
 *
 * Calendar time is cal BP: larger values are older, and a_start > b_end.
 */
functions {
#include susq_functions.stan
}

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> P;                                  // equals S here
  array[N] int<lower=1, upper=S> site;
  array[S] int<lower=1, upper=P> phase_start;
  array[S] int<lower=1> n_phase;
  array[N] int<lower=1, upper=3> class_id;

  int<lower=1> V;
  array[N] int<lower=1, upper=V> var_start;
  array[N] int<lower=1> n_var;
  array[V] int<lower=0> kappa_idx;
  array[V] int<lower=0, upper=1> is_outlier;

  int<lower=1> Ltot;
  vector[Ltot] L_flat;
  vector[Ltot] Phi_flat;
  array[V] int<lower=1> curve_pos;
  array[V] int<lower=1> curve_len;
  array[V] real curve_t0;
  array[V] real phi_total;
  real<lower=0> dt;

  int<lower=1> n_kappa;
  real cal_min;
  real cal_max;
  real mu_d_prior_mean;
  int<lower=0, upper=1> prior_only;
}

parameters {
  // 2S boundaries, increasing in cal BP: b_1 < a_1 < b_2 < a_2 < ...
  // so phase 1 is the youngest and phase S the oldest, and none overlap.
  // The calendar range is imposed by the declaration, which gives a proper
  // uniform prior over ordered boundaries inside it.
  ordered[2 * S] bound;
  real mu_d;
  real<lower=0> sigma_d;
  array[2] simplex[n_kappa] q;
  vector<lower=0, upper=1>[3] rho;
}

transformed parameters {
  vector[S] b_end;
  vector[S] a_start;
  vector[S] dur;
  for (s in 1:S) {
    b_end[s] = bound[2 * s - 1];
    a_start[s] = bound[2 * s];
    dur[s] = a_start[s] - b_end[s];
  }
}

model {
  // bound is uniform over its declared ordered range. The map from
  // (b_end, a_start) to (b_end, dur) is linear with unit Jacobian, so the
  // duration prior below needs no adjustment.
  dur ~ lognormal(mu_d, sigma_d);
  mu_d ~ normal(mu_d_prior_mean, 1);
  sigma_d ~ normal(0, 1);
  rho ~ beta(2, 18);
  for (c in 1:2) q[c] ~ dirichlet(rep_vector(1.0, n_kappa));

  if (prior_only == 0) {
    for (i in 1:N) {
      target += variant_logprob(a_start[site[i]], b_end[site[i]], i, class_id[i],
                                var_start, n_var, kappa_idx, is_outlier, rho, q,
                                L_flat, Phi_flat, curve_pos, curve_len,
                                curve_t0, phi_total, dt);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  array[S] real phase_start_calBP;
  array[S] real phase_end_calBP;
  array[S] real phase_duration;
  // Every phase holds its site's determinations by construction, so all are
  // occupation episodes. The name and shape match susq_phases.stan so the
  // same R metrics code reads either fit.
  array[S] int<lower=0, upper=1> active = rep_array(1, S);

  for (s in 1:S) {
    phase_start_calBP[s] = a_start[s];
    phase_end_calBP[s] = b_end[s];
    phase_duration[s] = dur[s];
  }
  for (i in 1:N) {
    log_lik[i] = variant_logprob(a_start[site[i]], b_end[site[i]], i, class_id[i],
                                 var_start, n_var, kappa_idx, is_outlier, rho, q,
                                 L_flat, Phi_flat, curve_pos, curve_len,
                                 curve_t0, phi_total, dt);
  }
}
