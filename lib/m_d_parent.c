#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern int m_d_parent(mpq_t angle, complex double *root_out, complex double *parent_out, complex double nucleus, int period, int maxsteps) {
  switch (m_d_shape(nucleus, period)) {
    case m_cardioid: {
      // find root directly
      complex double z = nucleus;
      complex double c = nucleus;
      m_d_interior(&z, &c, z, c, 1, period, maxsteps);
      *root_out = c;
      return 0; // no parent
    }
    case m_circle: {
      // trace internal ray to near to root
      complex double z = nucleus;
      complex double c = nucleus;
      for (int step = 0; step < maxsteps - 1; ++step) {
        m_d_interior(&z, &c, z, c, (step + 0.5) / maxsteps, period, maxsteps);
      }
      complex double root0 = c;
      m_d_interior(&z, &c, z, c, (maxsteps - 0.5) / maxsteps, period, maxsteps);
      complex double root1 = c;
      // find interior coordinate of a point just past the root into the parent
      complex double parent_guess = 2 * root1 - root0;
      c = parent_guess;
      z = 0;
      double mz2 = 1.0 / 0.0;
      for (int p = 1; p < period; ++p) {
        z = z * z + c;
        double z2 = cabs2(z);
        if (z2 < mz2) {
          mz2 = z2;
          if (period % p == 0) {
            complex double w = z;
            m_d_wucleus(&w, w, c, p, maxsteps);
            complex double dw = 1;
            for (int q = 0 ; q < p; ++q) {
              dw = 2 * w * dw;
              w = w * w + c;
            }
            if (cabs2(dw) < 1) {
              // interior to component of period p
              int den = period / p;
              int num = ((int) round(den * carg(w) / twopi) + den) % den;
              mpq_set_si(angle, num, den);
              mpq_canonicalize(angle);
              m_d_nucleus(&c, c, p, maxsteps);
              *parent_out = c;
              complex double interior = cexp(I * twopi * num / (double) den);
              m_d_interior(&w, &c, w, c, interior, p, maxsteps);
              *root_out = c;
              return p; // period of parent
            }
          }
        }
      }
    }
  }
  return -1; // fail
}
