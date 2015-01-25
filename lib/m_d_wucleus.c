#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_wucleus(complex double *z, complex double z_guess, complex double c, int period) {
  complex double zz = z_guess;
  complex double dzz = 1;
  for (int i = 0; i < period; ++i) {
    dzz = 2 * zz * dzz;
    zz = zz * zz + c;
  }
  if (cabs2(zz - z_guess) <= epsilon2) {
    *z = z_guess;
    return m_converged;
  }
  complex double z_new = z_guess - (zz - z_guess) / (dzz - 1);
  complex double d = z_new - zz;
  if (cabs2(d) <= epsilon2) {
    *z = z_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *z = z_new;
    return m_stepped;
  } else {
    *z = z_guess;
    return m_failed;
  }
}
