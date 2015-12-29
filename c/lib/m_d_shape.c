#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_shape m_d_shape(double _Complex nucleus, int period) {
  double _Complex z = nucleus;
  double _Complex dc = 1;
  double _Complex dz = 1;
  double _Complex dcdc = 0;
  double _Complex dcdz = 0;
  for (int i = 1; i < period; ++i) {
    dcdc = 2 * (z * dcdc + dc * dc);
    dcdz = 2 * (z * dcdz + dc * dz);
    dc = 2 * z * dc + 1;
    dz = 2 * z * dz;
    z = z * z + nucleus;
  }
  double _Complex e = - (dcdc / (2 * dc) + dcdz / dz) / (dc * dz);
  if (cabs2(e) < cabs2(e - 1)) {
    return m_cardioid;
  } else {
    return m_circle;
  }
}
