#include <mandelbrot-numerics.h>

extern m_newton m_r_wucleus_step_raw(mpc_t z_out, const mpc_t z_guess, const mpc_t c, int period, mpc_t z, mpc_t dz, mpc_t z_new, mpc_t d, mpfr_t d2, mpfr_t epsilon2) {
  // z = z_guess; dz = 1;
  mpc_set(z, z_guess, MPC_RNDNN);
  mpc_set_si(dz, 1, MPC_RNDNN);
  for (int i = 0; i < period; ++i) {
    // dz = 2 * z * dz;
    mpc_mul(dz, z, dz, MPC_RNDNN);
    mpc_mul_2si(dz, dz, 1, MPC_RNDNN);
    // z = z * z + c;
    mpc_sqr(z, z, MPC_RNDNN);
    mpc_add(z, z, c, MPC_RNDNN);
  }
  // z_new = z_guess - (z - z_guess) / (dz - 1);
  mpc_sub(z, z, z_guess, MPC_RNDNN);
  mpc_sub_ui(dz, dz, 1, MPC_RNDNN);
  mpc_div(z, z, dz, MPC_RNDNN);
  mpc_sub(z_new, z_guess, z, MPC_RNDNN);
  // d = z_new - z_guess;
  mpc_sub(d, z_new, z_guess, MPC_RNDNN);
  // d2 = norm(d);
  mpc_norm(d2, d, MPFR_RNDN);
  if (mpfr_lessequal_p(d2, epsilon2)) {
    mpc_set(z_out, z_new, MPC_RNDNN);
    return m_converged;
  }
  if (mpfr_number_p(mpc_realref(d)) && mpfr_number_p(mpc_imagref(d))) {
    mpc_set(z_out, z_new, MPC_RNDNN);
    return m_stepped;
  } else {
    mpc_set(z_out, z_guess, MPC_RNDNN);
    return m_failed;
  }
}

extern m_newton m_r_wucleus_step(mpc_t z_out, const mpc_t z_guess, const mpc_t c, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, z_guess);
  prec = precr > preci ? precr : preci;
  mpc_get_prec2(&precr, &preci, c);
  prec = precr > prec ? precr : prec;
  prec = preci > prec ? preci : prec;
  mpc_set_prec(z_out, prec); // FIXME might trash when z_out = z_guess
  // init
  mpc_t z, dz, z_new, d;
  mpfr_t d2, epsilon2;
  mpc_init2(z, prec);
  mpc_init2(dz, prec);
  mpc_init2(z_new, prec);
  mpc_init2(d, prec);
  mpfr_init2(d2, prec);
  mpfr_init2(epsilon2, prec);
  // epsilon
  mpfr_set_si(epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(epsilon2);
  mpfr_sub_si(epsilon2, epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(epsilon2, epsilon2, MPFR_RNDN);
  // step raw
  m_newton retval = m_r_wucleus_step_raw(z_out, z_guess, c, period, z, dz, z_new, d, d2, epsilon2);
  // cleanup
  mpc_clear(z);
  mpc_clear(dz);
  mpc_clear(z_new);
  mpc_clear(d);
  mpfr_clear(d2);
  mpfr_clear(epsilon2);
  return retval;
}

extern m_newton m_r_wucleus(mpc_t z_out, const mpc_t z_guess, const mpc_t c, int period, int maxsteps) {
  m_newton result = m_failed;
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, z_guess);
  prec = precr > preci ? precr : preci;
  mpc_get_prec2(&precr, &preci, c);
  prec = precr > prec ? precr : prec;
  prec = preci > prec ? preci : prec;
  // init
  mpc_t z0, z, dz, z_new, d;
  mpfr_t d2, epsilon2;
  mpc_init2(z0, prec);
  mpc_init2(z, prec);
  mpc_init2(dz, prec);
  mpc_init2(z_new, prec);
  mpc_init2(d, prec);
  mpfr_init2(d2, prec);
  mpfr_init2(epsilon2, prec);
  // epsilon
  mpfr_set_si(epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(epsilon2);
  mpfr_sub_si(epsilon2, epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(epsilon2, epsilon2, MPFR_RNDN);
  // z0 = z_guess;
  mpc_set(z0, z_guess, MPC_RNDNN);
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_r_wucleus_step_raw(z0, z0, c, period, z, dz, z_new, d, d2, epsilon2))) {
      break;
    }
  }
  // z_out = z0;
  mpc_set_prec(z_out, prec);
  mpc_set(z_out, z0, MPC_RNDNN);
  // cleanup
  mpc_clear(z0);
  mpc_clear(z);
  mpc_clear(dz);
  mpc_clear(z_new);
  mpc_clear(d);
  mpfr_clear(d2);
  mpfr_clear(epsilon2);
  return result;
}
