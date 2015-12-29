#include <mandelbrot-numerics.h>

extern double _Complex m_d_size(double _Complex nucleus, int period) {
  double _Complex l = 1;
  double _Complex b = 1;
  double _Complex z = 0;
  for (int i = 1; i < period; ++i) {
    z = z * z + nucleus;
    l = 2 * z * l;
    b = b + 1 / l;
  }
  return 1 / (b * l * l);
}
