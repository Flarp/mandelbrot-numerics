#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_nucleus_step(complex double *c_out, complex double c_guess, int period) {
  complex double z = 0;
  complex double dc = 0;
  for (int i = 0; i < period; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
  if (cabs2(dc) <= epsilon2) {
    *c_out = c_guess;
    return m_converged;
  }
  complex double c_new = c_guess - z / dc;
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

extern m_newton m_d_nucleus(complex double *c_out, complex double c_guess, int period, int maxsteps) {
  m_newton result = m_failed;
  complex double c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_nucleus_step(&c, c, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}
