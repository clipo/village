/**
 * Shared functions for the Susquehanna occupation models.
 *
 * Calendar time is cal BP throughout: larger values are older. A phase runs
 * from `a` (start, older, larger) to `b` (end, younger, smaller), a > b.
 */

/**
 * Cumulative integral of a variant's likelihood curve, evaluated at an
 * arbitrary calendar date.
 *
 * Each variant stores its likelihood L on a uniform grid of `curve_len[v]`
 * knots starting at cal BP `curve_t0[v]`, together with the exact integral
 * Phi of the piecewise-linear interpolant of L at each knot. Because L is
 * treated as piecewise linear, Phi is piecewise quadratic and therefore
 * continuously differentiable in t, which is what Hamiltonian Monte Carlo
 * needs. Outside the stored window the curve carries no meaningful mass, so
 * Phi clamps to 0 below and to phi_total above.
 */
real phi_at(real t, int v,
            vector L_flat, vector Phi_flat,
            array[] int curve_pos, array[] int curve_len,
            array[] real curve_t0, array[] real phi_total,
            real dt) {
  real lo = curve_t0[v];
  real hi = curve_t0[v] + (curve_len[v] - 1) * dt;
  if (t <= lo) return 0;
  if (t >= hi) return phi_total[v];

  // Locate the grid cell containing t. Stan cannot cast a parameter-dependent
  // real to an integer, so the index cannot be computed arithmetically even
  // though the grid is uniform. Branching on a parameter is permitted, so a
  // binary search finds the cell in about eleven comparisons. The selected
  // index is discrete but the expression below is continuously differentiable
  // across cell boundaries, so the gradient is correct and continuous.
  int g_lo = 1;
  int g_hi = curve_len[v];
  while (g_hi - g_lo > 1) {
    int mid = (g_lo + g_hi) %/% 2;
    if (lo + (mid - 1) * dt <= t) {
      g_lo = mid;
    } else {
      g_hi = mid;
    }
  }
  int g = g_lo;

  real h = t - (lo + (g - 1) * dt);
  int p = curve_pos[v] + g - 1;
  return Phi_flat[p] + L_flat[p] * h
         + (L_flat[p + 1] - L_flat[p]) * h * h / (2 * dt);
}

/**
 * Probability of a determination given a uniform occupation phase [b, a],
 * with the latent calendar date integrated out. `a` is the older boundary.
 * The 1/(a - b) factor is the uniform-phase normalisation, applied once per
 * determination, which is what makes phase duration identifiable.
 */
real phase_average(real a, real b, int v,
                   vector L_flat, vector Phi_flat,
                   array[] int curve_pos, array[] int curve_len,
                   array[] real curve_t0, array[] real phi_total,
                   real dt) {
  return (phi_at(a, v, L_flat, Phi_flat, curve_pos, curve_len,
                 curve_t0, phi_total, dt)
        - phi_at(b, v, L_flat, Phi_flat, curve_pos, curve_len,
                 curve_t0, phi_total, dt)) / (a - b);
}

/**
 * Log probability of determination `i` under one occupation phase, summing
 * over that determination's likelihood variants.
 *
 * Variants are the inbuilt-age scales plus one outlier variant. Their weights
 * are: rho[c] for the outlier; (1 - rho[c]) * q[c-1][kappa_idx] for an
 * inbuilt-age variant of wood (c = 2) or indeterminate (c = 3) material; and
 * (1 - rho[c]) for short-lived material (c = 1), which has a single variant
 * fixed at kappa = 0.
 *
 * The model block and the generated quantities block both need this, and two
 * copies could drift so that log_lik silently stopped matching the fitted
 * target, corrupting every LOO comparison. One implementation, called twice.
 */
real variant_logprob(real a, real b, int i, int c,
                     array[] int var_start, array[] int n_var,
                     array[] int kappa_idx, array[] int is_outlier,
                     vector rho, array[] vector q,
                     vector L_flat, vector Phi_flat,
                     array[] int curve_pos, array[] int curve_len,
                     array[] real curve_t0, array[] real phi_total,
                     real dt) {
  vector[n_var[i]] lp_var;
  for (m in 1:n_var[i]) {
    int v = var_start[i] + m - 1;
    // Every mixture log-weight needs the same protection as the likelihood
    // floor. A simplex component driven into a corner underflows to exactly
    // zero, log(0) is -inf, and its gradient is 1/0 = inf; log_sum_exp then
    // multiplies that infinity by a weight of essentially zero and yields NaN,
    // which diverges the transition. The inbuilt-age simplex q is the likely
    // offender here, since the data plainly prefer zero inbuilt age for most
    // material and push the other five components toward the corner.
    real lw;
    if (is_outlier[v] == 1) {
      lw = log(rho[c] + 1e-300);
    } else if (c == 1) {
      lw = log1m(rho[c]);
    } else {
      lw = log1m(rho[c]) + log(q[c - 1][kappa_idx[v]] + 1e-300);
    }
    real A = phase_average(a, b, v, L_flat, Phi_flat, curve_pos, curve_len,
                           curve_t0, phi_total, dt);
    // A is zero when the phase misses this variant's support entirely, which
    // is the common case once a site carries several phases. The floor must
    // not be tiny: the gradient of log(A + eps) is (dA/da) / (A + eps), so an
    // eps of 1e-300 lets the gradient overflow to infinity, and log_sum_exp
    // then multiplies that infinity by a mixture weight of essentially zero,
    // giving NaN and diverging the transition. Likelihoods here are normalised
    // so a phase average is of order 1e-2 to 1e-3, so 1e-12 is nine orders
    // below anything meaningful while keeping the gradient bounded.
    lp_var[m] = lw + log(A + 1e-12);
  }
  return log_sum_exp(lp_var);
}
