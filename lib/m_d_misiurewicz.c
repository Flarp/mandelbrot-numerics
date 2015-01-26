#include <mandelbrot-numerics.h>
#include "m_d_util.h"

/*
(pp+p) = (pp)
(pp+p) / (pp) - 1 = 0
  quotient rule:
    f(x) = \frac{g(x)}{h(x)}
   f'(x) = \frac{g'(x)h(x) - g(x)h'(x)}{[h(x)]^2}
((pp+p)' (pp) - (pp+p) (pp)') / (pp)^2
*/

extern m_newton m_d_misiurewicz_step(complex double *c_out, complex double c_guess, int preperiod, int period) {
  complex double z = c_guess;
  complex double dc = 1;
  complex double zp = 0;
  complex double dcp = 0;
  for (int i = 0; i < preperiod + period; ++i) {
    if (i == preperiod) {
      zp = z;
      dcp = dc;
    }
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
  complex double c_new = c_guess - (z / zp - 1) / ((dc * zp - z * dcp) / (zp * zp));
  complex double d = c_new - c_guess;
  if (cabs2(d) <= epsilon2) {
    *c_out = c_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *c_out = c_new;
    return m_stepped;
  } else {
    *c_out = c_guess;
    return m_failed;
  }
}

extern m_newton m_d_misiurewicz(complex double *c_out, complex double c_guess, int preperiod, int period, int maxsteps) {
  m_newton result = m_failed;
  complex double c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_misiurewicz_step(&c, c, preperiod, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}
