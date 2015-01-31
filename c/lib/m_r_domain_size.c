#include <mandelbrot-numerics.h>

extern void m_r_domain_size(mpfr_t size, const mpc_t nucleus, int period) {
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, nucleus);
  prec = precr > preci ? precr : preci;
  // init
#define CVARS z, dc
#define RVARS zq2, zp2
  mpc_t CVARS;
  mpfr_t RVARS;
  mpfr_inits2(prec, RVARS, (mpfr_ptr) 0);
  mpc_init2(z, prec);
  mpc_init2(dc, prec);
  // z = nucleus; dc = 1; zq2 = norm(z);
  mpc_set(z, nucleus, MPC_RNDNN);
  mpc_set_si(dc, 1, MPC_RNDNN);
  mpc_norm(zq2, z, MPFR_RNDN);
  for (int q = 2; q <= period; ++q) {
    // dc = 2 * z * dc + 1;
    mpc_mul(dc, z, dc, MPC_RNDNN);
    mpc_mul_2si(dc, dc, 1, MPC_RNDNN);
    mpc_add_ui(dc, dc, 1, MPC_RNDNN);
    // z = z * z + nucleus;
    mpc_sqr(z, z, MPC_RNDNN);
    mpc_add(z, z, nucleus, MPC_RNDNN);
    // zp2 = norm(z);
    mpc_norm(zp2, z, MPFR_RNDN);
    if (q < period && mpfr_less_p(zp2, zq2)) {
      mpfr_set(zq2, zp2, MPFR_RNDN);
    }
  }
  // size = sqrt(zq2 / norm(dc));
  mpfr_set_prec(size, prec);
  mpc_norm(zp2, dc, MPFR_RNDN);
  mpfr_div(size, zq2, zp2, MPFR_RNDN);
  mpfr_sqrt(size, size, MPFR_RNDN);
  // cleanup
  mpc_clear(z);
  mpc_clear(dc);
  mpfr_clears(RVARS, (mpfr_ptr) 0);
#undef CVARS
#undef RVARS
}
