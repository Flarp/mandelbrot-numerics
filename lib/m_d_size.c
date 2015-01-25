#include <mandelbrot-numerics.h>

extern complex double m_d_size(complex double nucleus, int period) {
  complex double l = 1;
  complex double b = 1;
  complex double z = 0;
  for (int i = 1; i < period; ++i) {
    z = z * z + nucleus;
    l = 2 * z * l;
    b = b + 1 / l;
  }
  return 1 / (b * l * l);
}
