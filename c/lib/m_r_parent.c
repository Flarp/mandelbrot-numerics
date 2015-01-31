#include <mandelbrot-numerics.h>

extern int m_r_parent(mpq_t angle_out, mpc_t root_out, mpc_t parent_out, const mpc_t nucleus, int period, int maxsteps) {
  int retval = -1;
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, nucleus);
  prec = precr > preci ? precr : preci;
  switch (m_r_shape(nucleus, period)) {
    case m_cardioid: { // find root directly
      // init
      mpc_t z, c, i;
      mpc_init2(z, prec);
      mpc_init2(c, prec);
      mpc_init2(i, prec);
      // z = nucleus; c = nucleus; i = 1;
      mpc_set(z, nucleus, MPC_RNDNN);
      mpc_set(c, nucleus, MPC_RNDNN);
      mpc_set_si(i, 1, MPC_RNDNN);
      m_r_interior(z, c, z, c, i, period, maxsteps);
      // root_out = c;
      mpc_set_prec(root_out, prec);
      mpc_set(root_out, c, MPC_RNDNN);
      // cleanup
      mpc_clear(z);
      mpc_clear(c);
      mpc_clear(i);
      retval = 0; // no parent
      break;
    }
    case m_circle: { // trace internal ray to near to root
      // init
      mpc_t z, c, i, root0, root1, parent_guess, w0, w, dw, interior;
      mpfr_t mz2, z2, d, pi;
      mpc_init2(z, prec);
      mpc_init2(c, prec);
      mpc_init2(i, prec);
      mpc_init2(root0, prec);
      mpc_init2(root1, prec);
      mpc_init2(parent_guess, prec);
      mpc_init2(w0, prec);
      mpc_init2(w, prec);
      mpc_init2(dw, prec);
      mpc_init2(interior, prec);
      mpfr_init2(mz2, prec);
      mpfr_init2(z2, prec);
      mpfr_init2(d, prec);
      mpfr_init2(pi, prec);
      mpfr_const_pi(pi, MPFR_RNDN);
      // z = nucleus; c = nucleus;
      mpc_set(z, nucleus, MPC_RNDNN);
      mpc_set(c, nucleus, MPC_RNDNN);
      for (int step = 0; step < maxsteps - 1; ++step) {
        // i = (step + 0.5) / maxsteps;
        mpc_set_d(i, (step + 0.5) / maxsteps, MPC_RNDNN);
        m_r_interior(z, c, z, c, i, period, maxsteps);
      }
      // root0 = c;
      mpc_set(root0, c, MPC_RNDNN);
      // i = (maxsteps - 0.5) / maxsteps;
      mpc_set_d(i, (maxsteps - 0.5) / maxsteps, MPC_RNDNN);
      m_r_interior(z, c, z, c, i, period, maxsteps);
      // root1 = c;
      mpc_set(root1, c, MPC_RNDNN);
      // find interior coordinate of a point just past the root into the parent
      // parent_guess = 2 * root1 - root0;
      mpc_mul_2si(root1, root1, 1, MPC_RNDNN);
      mpc_sub(parent_guess, root1, root0, MPC_RNDNN);
      // c = parent_guess; z = 0; mz2 = 1.0 / 0.0;
      mpc_set(c, parent_guess, MPC_RNDNN);
      mpc_set_si(z, 0, MPC_RNDNN);
      mpfr_set_d(mz2, 1.0 / 0.0, MPFR_RNDN);
      for (int p = 1; p < period; ++p) {
        // z = z * z + c;
        mpc_sqr(z, z, MPC_RNDNN);
        mpc_add(z, z, c, MPC_RNDNN);
        // z2 = norm(z);
        mpc_norm(z2, z, MPFR_RNDN);
        if (mpfr_less_p(z2, mz2)) {
          // mz2 = z2;
          mpfr_set(mz2, z2, MPFR_RNDN);
          if (period % p == 0) {
            m_r_wucleus(w0, z, c, p, maxsteps);
            // w = w0; dw = 1;
            mpc_set(w, w0, MPC_RNDNN);
            mpc_set_si(dw, 1, MPC_RNDNN);
            for (int q = 0 ; q < p; ++q) {
              // dw = 2 * w * dw;
              mpc_mul(dw, w, dw, MPC_RNDNN);
              mpc_mul_2si(dw, dw, 1, MPC_RNDNN);
              // w = w * w + c;
              mpc_sqr(w, w, MPC_RNDNN);
              mpc_add(w, w, c, MPC_RNDNN);
            }
            // d = norm(dw) - 1;
            mpc_norm(d, dw, MPFR_RNDN);
            mpfr_sub_ui(d, d, 1, MPFR_RNDN);
            if (mpfr_sgn(d) <= 0) {
              // interior to component of period p
              int den = period / p;
              // num = ((int) round(den * carg(w0) / twopi) + den) % den;
              mpc_arg(d, w0, MPFR_RNDN);
              mpfr_div(d, d, pi, MPFR_RNDN);
              mpfr_div_2ui(d, d, 1, MPFR_RNDN);
              mpfr_mul_si(d, d, den, MPFR_RNDN);
              mpfr_round(d, d);
              int num = (mpfr_get_si(d, MPFR_RNDN) + den) % den;
              mpq_set_si(angle_out, num, den);
              mpq_canonicalize(angle_out);
              m_r_nucleus(c, c, p, maxsteps);
              // parent_out = c;
              mpc_set(parent_out, c, MPC_RNDNN);
              // interior = cexp(I * twopi * num / (double) den);
              mpfr_mul_d(d, pi, 2 * mpq_get_d(angle_out), MPFR_RNDN);
              mpfr_sin_cos(mpc_imagref(interior), mpc_realref(interior), d, MPFR_RNDN);
              m_r_interior(w, c, w0, c, interior, p, maxsteps);
              // root_out = c;
              mpc_set(root_out, c, MPC_RNDNN);
              retval = p; // period of parent
              break;
            }
          }
        }
      }
      // cleanup
      mpc_clear(z);
      mpc_clear(c);
      mpc_clear(i);
      mpc_clear(root0);
      mpc_clear(root1);
      mpc_clear(parent_guess);
      mpc_clear(w0);
      mpc_clear(w);
      mpc_clear(dw);
      mpc_clear(interior);
      mpfr_clear(mz2);
      mpfr_clear(z2);
      mpfr_clear(d);
      mpfr_clear(pi);
      break;
    }
  }
  return retval; // fail
}
