#include <mandelbrot-numerics.h>

extern m_newton m_r_interior_step_raw(mpc_t z_out, mpc_t c_out, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period, mpc_t z_new, mpc_t c_new, mpc_t c, mpc_t z, mpc_t dz, mpc_t dc, mpc_t dzdz, mpc_t dcdz, mpc_t det, mpc_t d, mpc_t dz1, mpfr_t d2z, mpfr_t d2c, mpfr_t epsilon2) {
  // c = c_guess; z = z_guess; dz = 1; dc = 0; dzdz = 0; dcdz = 0;
  mpc_set(c, c_guess, MPC_RNDNN);
  mpc_set(z, z_guess, MPC_RNDNN);
  mpc_set_si(dz, 1, MPC_RNDNN);
  mpc_set_si(dc, 0, MPC_RNDNN);
  mpc_set_si(dzdz, 0, MPC_RNDNN);
  mpc_set_si(dcdz, 0, MPC_RNDNN);
  for (int p = 0; p < period; ++p) {
    // dcdz = 2 * (z * dcdz + dc * dz);
    mpc_mul(dcdz, z, dcdz, MPC_RNDNN);
    mpc_mul(d, dc, dz, MPC_RNDNN);
    mpc_add(dcdz, dcdz, d, MPC_RNDNN);
    mpc_mul_2si(dcdz, dcdz, 1, MPC_RNDNN);
    // dzdz = 2 * (z * dzdz + dz * dz);
    mpc_mul(dzdz, z, dzdz, MPC_RNDNN);
    mpc_sqr(d, dz, MPC_RNDNN);
    mpc_add(dzdz, dzdz, d, MPC_RNDNN);
    mpc_mul_2si(dzdz, dzdz, 1, MPC_RNDNN);
    // dc = 2 * z * dc + 1;
    mpc_mul(dc, z, dc, MPC_RNDNN);
    mpc_mul_2si(dc, dc, 1, MPC_RNDNN);
    mpc_add_ui(dc, dc, 1, MPC_RNDNN);
    // dz = 2 * z * dz;
    mpc_mul(dz, z, dz, MPC_RNDNN);
    mpc_mul_2si(dz, dz, 1, MPC_RNDNN);
    // z = z * z + c;
    mpc_sqr(z, z, MPC_RNDNN);
    mpc_add(z, z, c, MPC_RNDNN);
  }
  // det = (dz - 1) * dcdz - dc * dzdz;
  mpc_sub_ui(dz1, dz, 1, MPC_RNDNN);
  mpc_mul(det, dz1, dcdz, MPC_RNDNN);
  mpc_mul(d, dc, dzdz, MPC_RNDNN);
  mpc_sub(det, det, d, MPC_RNDNN);
  // z_new = z_guess - (dcdz * (z - z_guess) - dc * (dz - interior)) / det;
  mpc_sub(z, z, z_guess, MPC_RNDNN);
  mpc_sub(dz, dz, interior, MPC_RNDNN);
  mpc_mul(dcdz, dcdz, z, MPC_RNDNN);
  mpc_mul(dc, dc, dz, MPC_RNDNN);
  mpc_sub(dcdz, dcdz, dc, MPC_RNDNN);
  mpc_div(dcdz, dcdz, det, MPC_RNDNN);
  mpc_sub(z_new, z_guess, dcdz, MPC_RNDNN);
  // c_new = c_guess - ((dz - 1) * (dz - interior) - dzdz * (z - z_guess)) / det;
  mpc_mul(dz1, dz1, dz, MPC_RNDNN);
  mpc_mul(dzdz, dzdz, z, MPC_RNDNN);
  mpc_sub(dzdz, dz1, dzdz, MPC_RNDNN);
  mpc_div(dzdz, dzdz, det, MPC_RNDNN);
  mpc_sub(c_new, c_guess, dzdz, MPC_RNDNN);
  if (mpfr_number_p(mpc_realref(z_new)) && mpfr_number_p(mpc_imagref(z_new)) &&
      mpfr_number_p(mpc_realref(c_new)) && mpfr_number_p(mpc_imagref(c_new))) {
    // z_out = z_new; d2z = norm(z_new - z_guess);
    mpc_sub(d, z_new, z_guess, MPC_RNDNN); // in this order in case z_out = z_guess
    mpc_set(z_out, z_new, MPC_RNDNN);
    mpc_norm(d2z, d, MPFR_RNDN);
    // c_out = c_new; d2c = norm(c_new - c_guess);
    mpc_sub(d, c_new, c_guess, MPC_RNDNN); // in this order in case c_out = c_guess
    mpc_set(c_out, c_new, MPC_RNDNN);
    mpc_norm(d2c, d, MPFR_RNDN);
    if (mpfr_lessequal_p(d2z, epsilon2) && mpfr_lessequal_p(d2c, epsilon2)) {
      return m_converged;
    } else {
      return m_stepped;
    }
  } else {
    // z_out = z_guess; c_out = c_guess;
    mpc_set(z_out, z_guess, MPC_RNDNN);
    mpc_set(c_out, c_guess, MPC_RNDNN);
    return m_failed;
  }
}

extern m_newton m_r_interior_step(mpc_t z_out, mpc_t c_out, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period) {
  m_newton result = m_failed;
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, z_guess);
  prec = precr > preci ? precr : preci;
  mpc_get_prec2(&precr, &preci, c_guess);
  prec = precr > prec ? precr : prec;
  prec = preci > prec ? preci : prec;
  mpc_set_prec(z_out, prec); // FIXME might trash when z_out = z_guess
  mpc_set_prec(c_out, prec); // FIXME might trash when c_out = c_guess
  // init
  mpc_t z_new, c_new, c, z, dz, dc, dzdz, dcdz, det, d, dz1;
  mpfr_t d2z, d2c, epsilon2;
  mpc_init2(z_new, prec);
  mpc_init2(c_new, prec);
  mpc_init2(c, prec);
  mpc_init2(z, prec);
  mpc_init2(dc, prec);
  mpc_init2(dz, prec);
  mpc_init2(dcdz, prec);
  mpc_init2(dzdz, prec);
  mpc_init2(det, prec);
  mpc_init2(d, prec);
  mpc_init2(dz1, prec);
  mpfr_init2(d2z, prec);
  mpfr_init2(d2c, prec);
  mpfr_init2(epsilon2, prec);
  // epsilon
  mpfr_set_si(epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(epsilon2);
  mpfr_sub_si(epsilon2, epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(epsilon2, epsilon2, MPFR_RNDN);
  // step
  result = m_r_interior_step_raw(z_out, c_out, z_guess, c_guess, interior, period, z_new, c_new, c, z, dz, dc, dzdz, dcdz, det, d, dz1, d2z, d2c, epsilon2);
  // cleanup
  mpc_clear(z_new);
  mpc_clear(c_new);
  mpc_clear(c);
  mpc_clear(z);
  mpc_clear(dc);
  mpc_clear(dz);
  mpc_clear(dcdz);
  mpc_clear(dzdz);
  mpc_clear(det);
  mpc_clear(d);
  mpc_clear(dz1);
  mpfr_clear(d2z);
  mpfr_clear(d2c);
  mpfr_clear(epsilon2);
  return result;
}

extern m_newton m_r_interior(mpc_t z_out, mpc_t c_out, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period, int maxsteps) {
  m_newton result = m_failed;
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, z_guess);
  prec = precr > preci ? precr : preci;
  mpc_get_prec2(&precr, &preci, c_guess);
  prec = precr > prec ? precr : prec;
  prec = preci > prec ? preci : prec;
  // init
  mpc_t z, c, z_new, c_new, cc, zz, dz, dc, dzdz, dcdz, det, d, dz1;
  mpfr_t d2z, d2c, epsilon2;
  mpc_init2(z, prec);
  mpc_init2(c, prec);
  mpc_init2(z_new, prec);
  mpc_init2(c_new, prec);
  mpc_init2(cc, prec);
  mpc_init2(zz, prec);
  mpc_init2(dc, prec);
  mpc_init2(dz, prec);
  mpc_init2(dcdz, prec);
  mpc_init2(dzdz, prec);
  mpc_init2(det, prec);
  mpc_init2(d, prec);
  mpc_init2(dz1, prec);
  mpfr_init2(d2z, prec);
  mpfr_init2(d2c, prec);
  mpfr_init2(epsilon2, prec);
  // epsilon
  mpfr_set_si(epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(epsilon2);
  mpfr_sub_si(epsilon2, epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(epsilon2, epsilon2, MPFR_RNDN);
  // z = z_guess; c = c_guess;
  mpc_set(z, z_guess, MPC_RNDNN);
  mpc_set(c, c_guess, MPC_RNDNN);
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_r_interior_step_raw(z, c, z, c, interior, period, z_new, c_new, cc, zz, dz, dc, dzdz, dcdz, det, d, dz1, d2z, d2c, epsilon2))) {
      break;
    }
  }
  // z_out = z; c_out = c;
  mpc_set_prec(z_out, prec);
  mpc_set(z_out, z, MPC_RNDNN);
  mpc_set_prec(c_out, prec);
  mpc_set(c_out, c, MPC_RNDNN);
  // cleanup
  mpc_clear(z);
  mpc_clear(c);
  mpc_clear(z_new);
  mpc_clear(c_new);
  mpc_clear(cc);
  mpc_clear(zz);
  mpc_clear(dc);
  mpc_clear(dz);
  mpc_clear(dcdz);
  mpc_clear(dzdz);
  mpc_clear(det);
  mpc_clear(d);
  mpc_clear(dz1);
  mpfr_clear(d2z);
  mpfr_clear(d2c);
  mpfr_clear(epsilon2);
  return result;
}
