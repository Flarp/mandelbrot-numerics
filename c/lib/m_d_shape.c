#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern double _Complex m_d_shape_estimate(double _Complex nucleus, int period) {
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
  double _Complex shape = - (dcdc / (2 * dc) + dcdz / dz) / (dc * dz);
  return shape;
}

extern m_shape m_d_shape_discriminant(double _Complex shape) {
  if (cabs(shape) < cabs(shape - 1)) {
    return m_cardioid;
  } else {
    return m_circle;
  }
}

extern m_shape m_d_shape(double _Complex nucleus, int period) {
  return m_d_shape_discriminant(m_d_shape_estimate(nucleus, period));
}
