#include <mandelbrot-numerics.h>

extern m_newton m_r_nucleus_step_raw(mpc_t c_out, const mpc_t c_guess, int period, mpc_t z, mpc_t dc, mpc_t c_new, mpc_t d, mpfr_t d2, mpfr_t epsilon2) {
  // z = 0; dc = 0;
  mpc_set_si(z, 0, MPC_RNDNN);
  mpc_set_si(dc, 0, MPC_RNDNN);
  for (int i = 0; i < period; ++i) {
    // dc = 2 * z * dc + 1;
    mpc_mul(dc, z, dc, MPC_RNDNN);
    mpc_mul_2si(dc, dc, 1, MPC_RNDNN);
    mpc_add_ui(dc, dc, 1, MPC_RNDNN);
    // z = z * z + c_guess;
    mpc_sqr(z, z, MPC_RNDNN);
    mpc_add(z, z, c_guess, MPC_RNDNN);
  }
  mpc_norm(d2, dc, MPFR_RNDN);
  if (mpfr_lessequal_p(d2, epsilon2)) {
    mpc_set(c_out, c_guess, MPC_RNDNN);
    return m_converged;
  }
  // c_new = c_guess - z / dc;
  mpc_div(z, z, dc, MPC_RNDNN);
  mpc_sub(c_new, c_guess, z, MPC_RNDNN);
  // d = c_new - c_guess;
  mpc_sub(d, c_new, c_guess, MPC_RNDNN);
  mpc_norm(d2, d, MPFR_RNDN);
  if (mpfr_lessequal_p(d2, epsilon2)) {
    mpc_set(c_out, c_new, MPC_RNDNN);
    return m_converged;
  }
  if (mpfr_number_p(mpc_realref(d)) && mpfr_number_p(mpc_imagref(d))) {
    mpc_set(c_out, c_new, MPC_RNDNN);
    return m_stepped;
  } else {
    mpc_set(c_out, c_guess, MPC_RNDNN);
    return m_failed;
  }
}

extern m_newton m_r_nucleus_step(mpc_t c_out, const mpc_t c_guess, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, c_guess);
  prec = precr > preci ? precr : preci;
  mpc_set_prec(c_out, prec);
  // init
  mpc_t z, dc, c_new, d;
  mpfr_t d2, epsilon2;
  mpc_init2(z, prec);
  mpc_init2(dc, prec);
  mpc_init2(c_new, prec);
  mpc_init2(d, prec);
  mpfr_init2(d2, prec);
  mpfr_init2(epsilon2, prec);
  // epsilon
  mpfr_set_si(epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(epsilon2);
  mpfr_sub_si(epsilon2, epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(epsilon2, epsilon2, MPFR_RNDN);
  // step raw
  m_newton retval = m_r_nucleus_step_raw(c_out, c_guess, period, z, dc, c_new, d, d2, epsilon2);
  // cleanup
  mpc_clear(z);
  mpc_clear(dc);
  mpc_clear(c_new);
  mpc_clear(d);
  mpfr_clear(d2);
  mpfr_clear(epsilon2);
  return retval;
}

extern m_newton m_r_nucleus(mpc_t c_out, const mpc_t c_guess, int period, int maxsteps) {
  m_newton result = m_failed;
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, c_guess);
  prec = precr > preci ? precr : preci;
  mpc_set_prec(c_out, prec);
  // init
  mpc_t c, z, dc, c_new, d;
  mpfr_t d2, epsilon2;
  mpc_init2(c, prec);
  mpc_init2(z, prec);
  mpc_init2(dc, prec);
  mpc_init2(c_new, prec);
  mpc_init2(d, prec);
  mpfr_init2(d2, prec);
  mpfr_init2(epsilon2, prec);
  // epsilon
  mpfr_set_si(epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(epsilon2);
  mpfr_sub_si(epsilon2, epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(epsilon2, epsilon2, MPFR_RNDN);
  // c = c_guess
  mpc_set(c, c_guess, MPC_RNDNN);
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_r_nucleus_step_raw(c, c, period, z, dc, c_new, d, d2, epsilon2))) {
      break;
    }
  }
  mpc_set(c_out, c, MPC_RNDNN);
  // cleanup
  mpc_clear(c);
  mpc_clear(z);
  mpc_clear(dc);
  mpc_clear(c_new);
  mpc_clear(d);
  mpfr_clear(d2);
  mpfr_clear(epsilon2);
  return result;
}
