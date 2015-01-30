#include <mandelbrot-numerics.h>

extern m_shape m_r_shape(const mpc_t nucleus, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, nucleus);
  prec = precr > preci ? precr : preci;
  // init
  mpc_t z, dc, dz, dcdc, dcdz, e, e1;
  mpfr_t e2, e12;
  mpc_init2(z, prec);
  mpc_init2(dc, prec);
  mpc_init2(dz, prec);
  mpc_init2(dcdc, prec);
  mpc_init2(dcdz, prec);
  mpc_init2(e, prec);
  mpc_init2(e1, prec);
  mpfr_init2(e2, prec);
  mpfr_init2(e12, prec);
  // z = nucleus; dc = 1; dz = 1; dcdc = 0; dcdz = 0;
  mpc_set(z, nucleus, MPC_RNDNN);
  mpc_set_si(dc, 1, MPC_RNDNN);
  mpc_set_si(dz, 1, MPC_RNDNN);
  mpc_set_si(dcdc, 0, MPC_RNDNN);
  mpc_set_si(dcdz, 0, MPC_RNDNN);
  for (int i = 1; i < period; ++i) {
    // dcdc = 2 * (z * dcdc + dc * dc);
    mpc_mul(dcdc, z, dcdc, MPC_RNDNN);
    mpc_sqr(e, dc, MPC_RNDNN);
    mpc_add(dcdc, dcdc, e, MPC_RNDNN);
    mpc_mul_2si(dcdc, dcdc, 1, MPC_RNDNN);
    // dcdz = 2 * (z * dcdz + dc * dz);
    mpc_mul(dcdz, z, dcdz, MPC_RNDNN);
    mpc_mul(e, dc, dz, MPC_RNDNN);
    mpc_add(dcdz, dcdz, e, MPC_RNDNN);
    mpc_mul_2si(dcdz, dcdz, 1, MPC_RNDNN);
    // dc = 2 * z * dc + 1;
    mpc_mul(dc, z, dc, MPC_RNDNN);
    mpc_mul_2si(dc, dc, 1, MPC_RNDNN);
    mpc_add_ui(dc, dc, 1, MPC_RNDNN);
    // dz = 2 * z * dz;
    mpc_mul(dz, z, dz, MPC_RNDNN);
    mpc_mul_2si(dz, dz, 1, MPC_RNDNN);
    // z = z * z + nucleus;
    mpc_sqr(z, z, MPC_RNDNN);
    mpc_add(z, z, nucleus, MPC_RNDNN);
  }
  // e = (dcdc / (2 * dc) + dcdz / dz) / (dc * dz);
  mpc_div(dcdc, dcdc, dc, MPC_RNDNN);
  mpc_div_2si(dcdc, dcdc, 1, MPC_RNDNN);
  mpc_div(dcdz, dcdz, dc, MPC_RNDNN);
  mpc_add(dcdc, dcdc, dcdz, MPC_RNDNN);
  mpc_mul(dc, dc, dz, MPC_RNDNN);
  mpc_div(e, dcdc, dc, MPC_RNDNN);
  // e1 = e + 1; e2 = norm(e); e12 = norm(e1);
  mpc_add_ui(e1, e, 1, MPC_RNDNN);
  mpc_norm(e2, e, MPFR_RNDN);
  mpc_norm(e12, e1, MPFR_RNDN);
  m_shape retval;
  if (mpfr_less_p(e2, e12)) {
    retval = m_cardioid;
  } else {
    retval = m_circle;
  }
  // cleanup
  mpc_clear(z);
  mpc_clear(dc);
  mpc_clear(dz);
  mpc_clear(dcdc);
  mpc_clear(dcdz);
  mpc_clear(e);
  mpc_clear(e1);
  mpfr_clear(e2);
  mpfr_clear(e12);
  return retval;
}
