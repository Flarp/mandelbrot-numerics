#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_wucleus_step(double _Complex *z, double _Complex z_guess, double _Complex c, int period) {
  double _Complex zz = z_guess;
  double _Complex dzz = 1;
  for (int i = 0; i < period; ++i) {
    dzz = 2 * zz * dzz;
    zz = zz * zz + c;
  }
  if (cabs2(zz - z_guess) <= epsilon2) {
    *z = z_guess;
    return m_converged;
  }
  double _Complex z_new = z_guess - (zz - z_guess) / (dzz - 1);
  double _Complex d = z_new - zz;
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

extern m_newton m_d_wucleus(double _Complex *z_out, double _Complex z_guess, double _Complex c, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex z = z_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_wucleus_step(&z, z, c, period))) {
      break;
    }
  }
  *z_out = z;
  return result;
}
