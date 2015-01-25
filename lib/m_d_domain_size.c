#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern double m_d_domain_size(complex double nucleus, int period) {
  complex double z = nucleus;
  complex double dc = 1;
  double zq2 = cabs2(z);
  for (int q = 2; q <= period; ++q) {
    dc = 2 * z * dc + 1;
    z = z * z + nucleus;
    double zp2 = cabs2(z);
    if (q < period && zp2 < zq2) {
      zq2 = zp2;
    }
  }
  return sqrt(zq2) / cabs(dc);
}
