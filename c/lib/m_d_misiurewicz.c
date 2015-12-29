#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern m_newton m_d_misiurewicz_naive_step(double _Complex *c_out, double _Complex c_guess, int preperiod, int period) {
  double _Complex z = 0;
  double _Complex dc = 0;
  double _Complex zp = 0;
  double _Complex dcp = 0;
  for (int i = 0; i < preperiod + 1 + period; ++i) {
    if (i == preperiod + 1) {
      zp = z;
      dcp = dc;
    }
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
  double _Complex f = z - zp;
  double _Complex df = dc - dcp;
  double _Complex c_new = c_guess - f / df;
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

extern m_newton m_d_misiurewicz_naive(double _Complex *c_out, double _Complex c_guess, int preperiod, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_misiurewicz_naive_step(&c, c, preperiod, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}

/*
(pp+p) = (pp)
(pp+p) / (pp) - 1 = 0
  quotient rule:
    f(x) = \frac{g(x)}{h(x)}
   f'(x) = \frac{g'(x)h(x) - g(x)h'(x)}{[h(x)]^2}
((pp+p)' (pp) - (pp+p) (pp)') / (pp)^2

(0 +p) /= ( 0)
(1 +p) /= ( 1)
...
(pp+p) =  (pp)

        (pp+p) - (pp)
------------------------------------ = 0
((0 +p) - ( 0))*((1 + p) - ( 1)) ...

  product rule:
    \frac{d}{dx} \left [ \prod_{i=1}^k f_i(x) \right ]
    = \sum_{i=1}^k \left(\frac{d}{dx} f_i(x) \prod_{j\ne i} f_j(x) \right)
    = \left( \prod_{i=1}^k f_i(x) \right) \left( \sum_{i=1}^k \frac{f'_i(x)}{f_i(x)} \right)


h = ((0 +p) - ( 0))*((1 + p) - ( 1)) ...
h' = h * sum (((i +p)'-(i)') / ((i +p) - (i))
(((pp+p)' - (pp)') * h - ((pp+p) - (pp)) * h') / h^2

*/

extern m_newton m_d_misiurewicz_step(double _Complex *c_out, double _Complex c_guess, int preperiod, int period) {
  // iteration
  double _Complex z = 0;
  double _Complex dc = 0;
  double _Complex zp[preperiod + 1 + period + 1];
  double _Complex dcp[preperiod + 1 + period + 1];
  for (int i = 0; i < preperiod + 1 + period; ++i) {
    zp[i] = z;
    dcp[i] = dc;
    dc = 2 * z * dc + 1;
    z = z * z + c_guess;
  }
  zp[preperiod + 1 + period] = z;
  dcp[preperiod + 1 + period] = dc;
  double _Complex h = 1;
  double _Complex dh = 0;
  // reject lower preperiods
  for (int i = 0; i < preperiod + 1; ++i) {
    double _Complex k = zp[i + period] - zp[i];
    h = h * k;
    dh = dh + (dcp[i + period] - dcp[i]) / k;
  }
/*
  // reject lower periods
  for (int i = 1; i < period; ++i) {
    for (int j = 0; j < period; ++j) {
      double _Complex k = zp[preperiod + 2 + ((j + i) % period)] - zp[preperiod + 2 + j];
      h = h * k;
      dh = dh + (dcp[preperiod + 2 + ((j + i) % period)] - dcp[preperiod + 2 + j]) / k;
    }
  }
*/
  // build function
  dh = dh * h;
  double _Complex g = zp[preperiod + 1 + period] - zp[preperiod + 1];
  double _Complex dg = dcp[preperiod + 1 + period] - dcp[preperiod + 1];
  double _Complex f = g / h;
  double _Complex df = (dg * h - g * dh) / (h * h);
  // newton step
  double _Complex c_new = c_guess - f / df;
  // check convergence
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

extern m_newton m_d_misiurewicz(double _Complex *c_out, double _Complex c_guess, int preperiod, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex c = c_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_misiurewicz_step(&c, c, preperiod, period))) {
      break;
    }
  }
  *c_out = c;
  return result;
}
