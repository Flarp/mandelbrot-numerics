#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_nucleus(complex double *c, complex double c_guess, int period) {
  complex double z = 0;
  complex double dc = 0;
  for (int i = 0; i < period; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
  if (cabs2(z) <= epsilon2) {
    *c = c_guess;
    return m_converged;
  }
  complex double c_new = c_guess - z / dc;
  complex double d = c_new - c_guess;
  if (cabs2(d) <= epsilon2) {
    *c = c_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *c = c_new;
    return m_stepped;
  } else {
    *c = c_guess;
    return m_failed;
  }
}
