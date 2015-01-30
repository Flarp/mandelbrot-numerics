#include <mandelbrot-numerics.h>

extern void m_r_size(mpc_t size, const mpc_t nucleus, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, nucleus);
  prec = precr > preci ? precr : preci;
  mpc_set_prec(size, prec);
  // init
  mpc_t z, l, l1, b;
  mpc_init2(z, prec);
  mpc_init2(l, prec);
  mpc_init2(l1, prec);
  mpc_init2(b, prec);
  // l = 1; b = 1; z = 0;
  mpc_set_si(l, 1, MPC_RNDNN);
  mpc_set_si(b, 1, MPC_RNDNN);
  mpc_set_si(z, 0, MPC_RNDNN);
  for (int i = 1; i < period; ++i) {
    // z = z * z + nucleus;
    mpc_sqr(z, z, MPC_RNDNN);
    mpc_add(z, z, nucleus, MPC_RNDNN);
    // l = 2 * z * l;
    mpc_mul(l, z, l, MPC_RNDNN);
    mpc_mul_2si(l, l, 1, MPC_RNDNN);
    // b = b + 1 / l;
    mpc_ui_div(l1, 1, l, MPC_RNDNN);
    mpc_add(b, b, l1, MPC_RNDNN);
  }
  // size = 1 / (b * l * l);
  mpc_sqr(l, l, MPC_RNDNN);
  mpc_mul(l1, b, l, MPC_RNDNN);
  mpc_ui_div(size, 1, l1, MPC_RNDNN);
  // cleanup
  mpc_clear(z);
  mpc_clear(l);
  mpc_clear(l1);
  mpc_clear(b);
}
