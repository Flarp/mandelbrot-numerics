#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_nucleus_step(double _Complex *c_out, double _Complex c_guess, int period) {
  double _Complex z = 0;
  double _Complex dc = 0;
  for (int i = 0; i < period; ++i) {
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
  if (cabs2(dc) <= epsilon2) {
    *c_out = c_guess;
    return m_converged;
  }
  double _Complex c_new = c_guess - z / dc;
  double _Complex d = c_new - c_guess;
  if (cabs2(d) <= epsilon2) {
    *c_out = c_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *c_out = c_new;
    return m_stepped;
  } else {
    *c_out = c_guess;
    return m_failed;
  }
}

extern m_newton m_d_nucleus(double _Complex *c_out, double _Complex c_guess, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_nucleus_step(&c, c, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}
