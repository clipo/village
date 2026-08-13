// Harness that exposes the interpolant to R for direct numerical comparison.
functions {
#include susq_functions.stan
}
data {
  int<lower=1> Ltot;
  vector[Ltot] L_flat;
  vector[Ltot] Phi_flat;
  int<lower=1> V;
  array[V] int curve_pos;
  array[V] int curve_len;
  array[V] real curve_t0;
  array[V] real phi_total;
  real dt;
  int<lower=1> M;
  array[M] real t_query;
  array[M] int v_query;
}
parameters { real dummy; }
model { dummy ~ std_normal(); }
generated quantities {
  array[M] real phi_out;
  for (m in 1:M)
    phi_out[m] = phi_at(t_query[m], v_query[m], L_flat, Phi_flat,
                        curve_pos, curve_len, curve_t0, phi_total, dt);
}
