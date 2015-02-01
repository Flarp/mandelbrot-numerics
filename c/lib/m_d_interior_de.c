#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern bool m_d_interior_de(double *de_out, complex double *dz_out, complex double z, complex double c, int p, int steps) {
  complex double z00 = 0;
  if (m_failed != m_d_wucleus(&z00, z, c, p, steps)) {
    complex double z0 = z00;
    complex double dz0 = 1;
    for (int j = 0; j < p; ++j) {
      dz0 = 2 * z0 * dz0;
      z0 = z0 * z0 + c;
    }
    if (cabs2(dz0) <= 1) {
      complex double z1 = z00;
      complex double dz1 = 1;
      complex double dzdz1 = 0;
      complex double dc1 = 0;
      complex double dcdz1 = 0;
      for (int j = 0; j < p; ++j) {
        dcdz1 = 2 * (z1 * dcdz1 + dz1 * dc1);
        dc1 = 2 * z1 * dc1 + 1;
        dzdz1 = 2 * (dz1 * dz1 + z1 * dzdz1);
        dz1 = 2 * z1 * dz1;
        z1 = z1 * z1 + c;
      }
      *de_out = (1 - cabs2(dz1)) / cabs(dcdz1 + dzdz1 * dc1 / (1 - dz1));
      *dz_out = dz1;
      return true;
    }
  }
  return false;
}
