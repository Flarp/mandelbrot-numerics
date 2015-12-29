#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_interior_step(double _Complex *z_out, double _Complex *c_out, double _Complex z_guess, double _Complex c_guess, double _Complex interior, int period) {
  double _Complex c = c_guess;
  double _Complex z = z_guess;
  double _Complex dz = 1;
  double _Complex dc = 0;
  double _Complex dzdz = 0;
  double _Complex dcdz = 0;
  for (int p = 0; p < period; ++p) {
    dcdz = 2 * (z * dcdz + dc * dz);
    dzdz = 2 * (z * dzdz + dz * dz);
    dc = 2 * z * dc + 1;
    dz = 2 * z * dz;
    z = z * z + c;
  }
  double _Complex det = (dz - 1) * dcdz - dc * dzdz;
  double _Complex z_new = z_guess - (dcdz * (z - z_guess) - dc * (dz - interior)) / det;
  double _Complex c_new = c_guess - ((dz - 1) * (dz - interior) - dzdz * (z - z_guess)) / det;
  if (cisfinite(z_new) && cisfinite(c_new)) {
    *z_out = z_new;
    *c_out = c_new;
    if (cabs2(z_new - z_guess) <= epsilon2 && cabs2(c_new - c_guess) <= epsilon2) {
      return m_converged;
    } else {
      return m_stepped;
    }
  } else {
    *z_out = z_guess;
    *c_out = c_guess;
    return m_failed;
  }
}

extern m_newton m_d_interior(double _Complex *z_out, double _Complex *c_out, double _Complex z_guess, double _Complex c_guess, double _Complex interior, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex z = z_guess;
  double _Complex c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_interior_step(&z, &c, z, c, interior, period))) {
      break;
    }
  }
  *z_out = z;
  *c_out = c;
  return result;
}
