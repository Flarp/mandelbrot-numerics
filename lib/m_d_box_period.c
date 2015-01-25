#include <mandelbrot-numerics.h>
#include "m_d_util.h"

static double cross(complex double a, complex double b) {
  return cimag(a) * creal(b) - creal(a) * cimag(b);
}

static bool crosses_positive_real_axis(complex double a, complex double b) {
  if (sgn(cimag(a)) != sgn(cimag(b))) {
    complex double d = b - a;
    int s = sgn(cimag(d));
    int t = sgn(cross(d, a));
    return s == t;
  }
  return false;
}

static bool surrounds_origin(complex double a, complex double b, complex double c, complex double d) {
  return odd
    ( crosses_positive_real_axis(a, b)
    + crosses_positive_real_axis(b, c)
    + crosses_positive_real_axis(c, d)
    + crosses_positive_real_axis(d, a)
    );
}

struct m_d_box_period {
  complex double c[4];
  complex double z[4];
  int p;
};

extern m_d_box_period *m_d_box_period_new(complex double center, double radius) {
  m_d_box_period *box = malloc(sizeof(*box));
  if (! box) {
    return 0;
  }
  box->z[0] = box->c[0] = center + ((-radius) + I * (-radius));
  box->z[1] = box->c[1] = center + (( radius) + I * (-radius));
  box->z[2] = box->c[2] = center + (( radius) + I * ( radius));
  box->z[3] = box->c[3] = center + ((-radius) + I * ( radius));
  box->p = 1;
  return box;
}

extern void m_d_box_period_delete(m_d_box_period *box) {
  if (box) {
    free(box);
  }
}

extern bool m_d_box_period_step(m_d_box_period *box) {
  if (! box) {
    return false;
  }
  bool ok = true;
  for (int i = 0; i < 4; ++i) {
    box->z[i] = box->z[i] * box->z[i] + box->c[i];
    ok = ok && cisfinite(box->z[i]);
  }
  box->p = box->p + 1;
  return ok;
}

extern bool m_d_box_period_have_period(const m_d_box_period *box) {
  if (! box) {
    return true;
  }
  return surrounds_origin(box->z[0], box->z[1], box->z[2], box->z[3]);
}

extern int m_d_box_period_get_period(const m_d_box_period *box) {
  if (! box) {
    return 0;
  }
  return box->p;
}

extern int m_d_box_period_do(complex double center, double radius, int maxperiod) {
  m_d_box_period *box = m_d_box_period_new(center, radius);
  if (! box) {
    return 0;
  }
  int period = 0;
  for (int i = 0; i < maxperiod; ++i) {
    if (m_d_box_period_have_period(box)) {
      period = m_d_box_period_get_period(box);
      break;
    }
    if (! m_d_box_period_step(box)) {
      break;
    }
  }
  m_d_box_period_delete(box);
  return period;
}
