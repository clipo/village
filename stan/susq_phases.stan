/**
 * Multi-phase occupation model for radiocarbon determinations grouped by site.
 *
 * Each site carries up to four candidate occupation phases. Determinations are
 * uniformly distributed within whichever phase they belong to, and phase
 * membership is marginalised as a Dirichlet mixture. The latent calendar date
 * of each determination is marginalised analytically using precomputed
 * cumulative calibration integrals, so no per-determination date parameter
 * appears. Inbuilt age and measurement outliers enter as a second mixture over
 * precomputed likelihood variants.
 *
 * Setting every n_phase to 1 gives the independent-phase model; collapsing all
 * determinations onto one site gives the single-phase model.
 *
 * Calendar time is cal BP: larger values are older, and a_start > b_end.
 */
functions {
#include susq_functions.stan
}

data {
  int<lower=1> N;                                  // determinations
  int<lower=1> S;                                  // sites
  int<lower=1> P;                                  // candidate phases in total
  array[N] int<lower=1, upper=S> site;
  array[S] int<lower=1, upper=P> phase_start;
  array[S] int<lower=1> n_phase;
  array[N] int<lower=1, upper=3> class_id;         // 1 short_lived, 2 wood, 3 indeterminate

  int<lower=1> V;                                  // likelihood variants in total
  array[N] int<lower=1, upper=V> var_start;
  array[N] int<lower=1> n_var;
  array[V] int<lower=0> kappa_idx;                 // 0 marks the outlier variant
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

  int<lower=0> n_site_J1;
  array[n_site_J1] int<lower=1, upper=S> sites_J1;
  int<lower=0> n_site_J2;
  array[n_site_J2] int<lower=1, upper=S> sites_J2;
  int<lower=0> n_site_J3;
  array[n_site_J3] int<lower=1, upper=S> sites_J3;
  int<lower=0> n_site_J4;
  array[n_site_J4] int<lower=1, upper=S> sites_J4;

  real mu_d_prior_mean;
  // Dirichlet concentration for the phase weights. The spec sets this to
  // alpha / J with alpha = 1, which is what makes the mixture overfitted so
  // that unsupported components empty out. Exposed as data so the sensitivity
  // of the sampler's geometry to it can be measured rather than assumed.
  // Concentration of the Dirichlet on phase weights. Below 1 the density is
  // unbounded at the simplex boundary, so the sampler is pulled toward weights
  // of exactly zero, where log(pi) underflows to negative infinity. That
  // sparsity pressure was needed only while spare components had to be emptied
  // out; with phase counts set from the data, every component is occupied and
  // a concentration of 1 or more is both sufficient and well behaved.
  // 1 pools phase durations across sites through mu_d and sigma_d; 0 gives
  // each phase an independent lognormal(mu_d_prior_mean, fixed_sigma_d) prior.
  // Pooling helps thinly dated phases but couples 73 durations of very
  // different identifiability to two hyperparameters, which is a funnel.
  int<lower=0, upper=1> hier_duration;
  real<lower=0> fixed_sigma_d;
  real<lower=0> dirichlet_conc_J2;
  real<lower=0> dirichlet_conc_J3;
  real<lower=0> dirichlet_conc_J4;
  int<lower=0, upper=1> prior_only;                // 1 runs the prior predictive
}

transformed data {
  // Map each site to its slot in the J=2 or J=4 simplex arrays. Sites with a
  // single phase need no simplex.
  array[S] int slot1 = rep_array(0, S);
  array[S] int slot2 = rep_array(0, S);
  array[S] int slot3 = rep_array(0, S);
  array[S] int slot4 = rep_array(0, S);
  for (m in 1:n_site_J1) slot1[sites_J1[m]] = m;
  for (m in 1:n_site_J2) slot2[sites_J2[m]] = m;
  for (m in 1:n_site_J3) slot3[sites_J3[m]] = m;
  for (m in 1:n_site_J4) slot4[sites_J4[m]] = m;
  real cal_range = cal_max - cal_min;
}

parameters {
  // Phase midpoints, ordered within each site and confined to the calendar
  // range. A site's J midpoints are the first J cumulative sums of a
  // (J+1)-simplex scaled across the range, so they are automatically ordered
  // and automatically inside it. A Dirichlet(1) on the simplex makes them
  // uniform over the grid subject to that ordering, which is what the spec
  // asks for. The earlier construction bounded only the first midpoint and
  // reached the rest through unbounded gaps, which left the parameters of
  // data-free mixture components on flat unbounded directions and produced
  // divergences in 70% of transitions.
  vector<lower=cal_min, upper=cal_max>[n_site_J1] mid1;
  array[n_site_J2] simplex[3] w2;
  array[n_site_J3] simplex[4] w3;
  array[n_site_J4] simplex[5] w4;

  // Durations, non-centred on a hierarchical lognormal.
  vector[P] log_dur_raw;
  array[hier_duration] real mu_d_param;
  array[hier_duration] real<lower=0> sigma_d_param;

  // Phase weights.
  array[n_site_J2] simplex[2] pi2;
  array[n_site_J3] simplex[3] pi3;
  array[n_site_J4] simplex[4] pi4;

  // Inbuilt-age distributions for wood and indeterminate material.
  array[2] simplex[n_kappa] q;

  // Outlier rate per material class.
  vector<lower=0, upper=1>[3] rho;
}

transformed parameters {
  vector[P] mid;
  real mu_d = hier_duration ? mu_d_param[1] : mu_d_prior_mean;
  real sigma_d = hier_duration ? sigma_d_param[1] : fixed_sigma_d;
  // Floor the duration. exp() underflows to exactly zero once its argument
  // drops below about -745, which the prior can reach while warming up, and a
  // zero-width phase makes the phase average 0/0. A thousandth of a year is
  // far below any archaeologically meaningful width.
  vector[P] dur = fmax(exp(mu_d + sigma_d * log_dur_raw), 1e-3);
  vector[P] a_start;
  vector[P] b_end;
  for (s in 1:S) {
    int p0 = phase_start[s];
    if (n_phase[s] == 1) {
      mid[p0] = mid1[slot1[s]];
    } else if (n_phase[s] == 2) {
      real c = 0;
      for (k in 1:2) { c += w2[slot2[s]][k]; mid[p0 + k - 1] = cal_min + cal_range * c; }
    } else if (n_phase[s] == 3) {
      real c = 0;
      for (k in 1:3) { c += w3[slot3[s]][k]; mid[p0 + k - 1] = cal_min + cal_range * c; }
    } else {
      real c = 0;
      for (k in 1:4) { c += w4[slot4[s]][k]; mid[p0 + k - 1] = cal_min + cal_range * c; }
    }
  }
  a_start = mid + 0.5 * dur;
  b_end   = mid - 0.5 * dur;
}

model {
  // Dirichlet(1) on the spacing simplex gives midpoints uniform over the
  // calendar range subject to their within-site ordering.
  for (m in 1:n_site_J2) w2[m] ~ dirichlet(rep_vector(1.0, 3));
  for (m in 1:n_site_J3) w3[m] ~ dirichlet(rep_vector(1.0, 4));
  for (m in 1:n_site_J4) w4[m] ~ dirichlet(rep_vector(1.0, 5));
  log_dur_raw ~ std_normal();
  if (hier_duration) {
    mu_d_param[1] ~ normal(mu_d_prior_mean, 1);
    sigma_d_param[1] ~ normal(0, 1);  // half-normal through the lower bound
  }
  rho ~ beta(2, 18);
  // alpha / J, with alpha = 1: an overfitted mixture, so phases unsupported by
  // data are driven toward zero weight rather than splitting the data.
  for (m in 1:n_site_J2) pi2[m] ~ dirichlet(rep_vector(dirichlet_conc_J2, 2));
  for (m in 1:n_site_J3) pi3[m] ~ dirichlet(rep_vector(dirichlet_conc_J3, 3));
  for (m in 1:n_site_J4) pi4[m] ~ dirichlet(rep_vector(dirichlet_conc_J4, 4));
  for (c in 1:2) q[c] ~ dirichlet(rep_vector(1.0, n_kappa));

  // The model block cannot return early, so the likelihood is guarded.
  if (prior_only == 0) {
    for (i in 1:N) {
      int s = site[i];
      int p0 = phase_start[s];
      int J = n_phase[s];
      vector[J] lp_phase;
      for (j in 1:J) {
        int p = p0 + j - 1;
        real lpi;
        if (J == 1) lpi = 0;
        else if (J == 2) lpi = log(pi2[slot2[s]][j] + 1e-300);
        else if (J == 3) lpi = log(pi3[slot3[s]][j] + 1e-300);
        else lpi = log(pi4[slot4[s]][j] + 1e-300);
        lp_phase[j] = lpi
          + variant_logprob(a_start[p], b_end[p], i, class_id[i],
                            var_start, n_var, kappa_idx, is_outlier, rho, q,
                            L_flat, Phi_flat, curve_pos, curve_len,
                            curve_t0, phi_total, dt);
      }
      target += log_sum_exp(lp_phase);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int z;                 // phase allocation, drawn from its conditional
  array[P] int<lower=0, upper=1> active = rep_array(0, P);
  array[P] real phase_start_calBP;
  array[P] real phase_end_calBP;
  array[P] real phase_duration;
  array[P] int n_allocated = rep_array(0, P);

  for (p in 1:P) {
    phase_start_calBP[p] = a_start[p];
    phase_end_calBP[p] = b_end[p];
    phase_duration[p] = dur[p];
  }

  for (i in 1:N) {
    int s = site[i];
    int p0 = phase_start[s];
    int J = n_phase[s];
    vector[J] lp_phase;
    for (j in 1:J) {
      int p = p0 + j - 1;
      real lpi;
      if (J == 1) lpi = 0;
      else if (J == 2) lpi = log(pi2[slot2[s]][j] + 1e-300);
      else if (J == 3) lpi = log(pi3[slot3[s]][j] + 1e-300);
      else lpi = log(pi4[slot4[s]][j] + 1e-300);
      lp_phase[j] = lpi
        + variant_logprob(a_start[p], b_end[p], i, class_id[i],
                          var_start, n_var, kappa_idx, is_outlier, rho, q,
                          L_flat, Phi_flat, curve_pos, curve_len,
                          curve_t0, phi_total, dt);
    }
    log_lik[i] = log_sum_exp(lp_phase);

    // Recover the discrete allocation the likelihood marginalised away.
    z[i] = p0 + categorical_logit_rng(lp_phase) - 1;
    n_allocated[z[i]] += 1;
  }

  // A phase is an occupation episode when at least one determination is
  // allocated to it. An episode exists when there is dated evidence for it,
  // which needs no weight threshold.
  for (p in 1:P) active[p] = n_allocated[p] > 0 ? 1 : 0;
}
