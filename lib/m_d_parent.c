#include <mandelbrot-numerics.h>
#include "m_d_util.h"

extern int m_d_parent(mpq_t angle, complex double *root_out, complex double *parent_out, complex double nucleus, int period, int maxsteps) {
  switch (m_d_shape(nucleus, period)) {
    case m_cardioid: {
      complex double z = nucleus;
      complex double c = nucleus;
      for (int i = 0; i < maxsteps; ++i) {
        if (m_stepped != m_d_interior(&z, &c, z, c, 1, period)) {
          break;
        }
      }
      *root_out = c;
      return 0;
    }
    case m_circle: {
      complex double z = nucleus;
      complex double c = nucleus;
      for (int step = 0; step < maxsteps - 1; ++step) {
        for (int i = 0; i < maxsteps; ++i) {
          if (m_stepped != m_d_interior(&z, &c, z, c, (step + 0.5) / maxsteps, period)) {
            break;
          }
        }
      }
      complex double root0 = c;
      for (int i = 0; i < maxsteps; ++i) {
        if (m_stepped != m_d_interior(&z, &c, z, c, (maxsteps - 0.5) / maxsteps, period)) {
          break;
        }
      }
      complex double root1 = c;
      complex double parent_guess = 2 * root1 - root0;
      // find interior coordinate
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
            for (int i = 0; i < maxsteps; ++i) {
              if (m_stepped != m_d_wucleus(&w, w, c, period)) {
                break;
              }
            }
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
              for (int i = 0; i < maxsteps; ++i) {
                if (m_stepped != m_d_nucleus(&c, c, period)) {
                  break;
                }
              }
              *parent_out = c;
              complex double interior = cexp(I * twopi * num / (double) den);
              for (int i = 0; i < maxsteps; ++i) {
                if (m_stepped != m_d_interior(&w, &c, w, c, interior, p)) {
                  break;
                }
              }
              *root_out = c;
              return p;
            }
          }
        }
      }
    }
  }
  return -1;
}
