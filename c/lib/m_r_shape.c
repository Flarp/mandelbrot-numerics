#include <mandelbrot-numerics.h>

extern void m_r_shape_estimate(mpc_t shape, const mpc_t nucleus, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, nucleus);
  prec = precr > preci ? precr : preci;
  // init
  mpc_t z, dc, dz, dcdc, dcdz, t;
  mpc_init2(z, prec);
  mpc_init2(dc, prec);
  mpc_init2(dz, prec);
  mpc_init2(dcdc, prec);
  mpc_init2(dcdz, prec);
  mpc_init2(t, prec);
  // z = nucleus; dc = 1; dz = 1; dcdc = 0; dcdz = 0;
  mpc_set(z, nucleus, MPC_RNDNN);
  mpc_set_si(dc, 1, MPC_RNDNN);
  mpc_set_si(dz, 1, MPC_RNDNN);
  mpc_set_si(dcdc, 0, MPC_RNDNN);
  mpc_set_si(dcdz, 0, MPC_RNDNN);
  for (int i = 1; i < period; ++i) {
    // dcdc = 2 * (z * dcdc + dc * dc);
    mpc_mul(dcdc, z, dcdc, MPC_RNDNN);
    mpc_sqr(t, dc, MPC_RNDNN);
    mpc_add(dcdc, dcdc, t, MPC_RNDNN);
    mpc_mul_2si(dcdc, dcdc, 1, MPC_RNDNN);
    // dcdz = 2 * (z * dcdz + dc * dz);
    mpc_mul(dcdz, z, dcdz, MPC_RNDNN);
    mpc_mul(t, dc, dz, MPC_RNDNN);
    mpc_add(dcdz, dcdz, t, MPC_RNDNN);
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
  // shape = -(dcdc / (2 * dc) + dcdz / dz) / (dc * dz);
  mpc_div(dcdc, dcdc, dc, MPC_RNDNN);
  mpc_div_2si(dcdc, dcdc, 1, MPC_RNDNN);
  mpc_div(dcdz, dcdz, dz, MPC_RNDNN);
  mpc_add(dcdc, dcdc, dcdz, MPC_RNDNN);
  mpc_mul(dc, dc, dz, MPC_RNDNN);
  mpc_div(shape, dcdc, dc, MPC_RNDNN);
  mpc_neg(shape, shape, MPC_RNDNN);
  // cleanup
  mpc_clear(z);
  mpc_clear(dc);
  mpc_clear(dz);
  mpc_clear(dcdc);
  mpc_clear(dcdz);
  mpc_clear(t);
}

extern m_shape m_r_shape_discriminant(const mpc_t shape) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, shape);
  prec = precr > preci ? precr : preci;
  // shape1 = shape - 1
  mpc_t one, shape1;
  mpc_init2(one, prec);
  mpc_init2(shape1, prec);
  mpc_set_si(one, 1, MPC_RNDNN);
  mpc_sub(shape1, shape, one, MPC_RNDNN);
  // e, e1
  mpfr_t e, e1;
  mpfr_init2(e, prec);
  mpfr_init2(e1, prec);
  // e = norm(shape); e1 = norm(shape1);
  mpc_norm(e, shape, MPFR_RNDN);
  mpc_norm(e1, shape1, MPFR_RNDN);
  m_shape retval;
  if (mpfr_less_p(e, e1)) {
    retval = m_cardioid;
  } else {
    retval = m_circle;
  }
  // cleanup
  mpc_clear(one);
  mpc_clear(shape1);
  mpfr_clear(e);
  mpfr_clear(e1);
  return retval;
}

extern m_shape m_r_shape(const mpc_t nucleus, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, nucleus);
  prec = precr > preci ? precr : preci;
  // shape
  mpc_t shape;
  mpc_init2(shape, prec);
  // retval = discriminant(estimate(nucleus, period))
  m_r_shape_estimate(shape, nucleus, period);
  m_shape retval = m_r_shape_discriminant(shape);
  // cleanup
  mpc_clear(shape);
  return retval;
}
