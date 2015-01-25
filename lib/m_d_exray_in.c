#include <mandelbrot-numerics.h>
#include "m_d_util.h"

struct m_d_exray_in {
  mpq_t angle;
  mpq_t one;
  int sharpness;
  double er;
  complex double c;
  int j;
  int k;
};

extern m_d_exray_in *m_d_exray_in_new(const mpq_t angle, int sharpness) {
  m_d_exray_in *ray = malloc(sizeof(*ray));
  mpq_init(ray->angle);
  mpq_set(ray->angle, angle);
  mpq_init(ray->one);
  mpq_set_ui(ray->one, 1, 1);
  ray->sharpness = sharpness;
  ray->er = 65536.0;
  double a = twopi * mpq_get_d(ray->angle);
  ray->c = ray->er * (cos(a) + I * sin(a));
  ray->k = 0;
  ray->j = 0;
  return ray;
}

extern void m_d_exray_in_delete(m_d_exray_in *ray) {
  if (! ray) {
    return;
  }
  mpq_clear(ray->angle);
  mpq_clear(ray->one);
  free(ray);
}

extern m_newton m_d_exray_in_step(m_d_exray_in *ray) {
  if (ray->j >= ray->sharpness) {
    mpq_mul_2exp(ray->angle, ray->angle, 1);
    if (mpq_cmp_ui(ray->angle, 1, 1) >= 0) {
      mpq_sub(ray->angle, ray->angle, ray->one);
    }
    ray->k = ray->k + 1;
    ray->j = 0;
  }
  double r = pow(ray->er, pow(0.5, (ray->j + 0.5) / ray->sharpness));
  double a = twopi * mpq_get_d(ray->angle);
  complex double target = r * (cos(a) + I * sin(a));
  complex double c = ray->c;
  for (int i = 0; i < 64; ++i) { // FIXME arbitrary limit
    complex double z = 0;
    complex double dc = 0;
    for (int p = 0; p <= ray->k; ++p) {
      dc = 2 * z * dc + 1;
      z = z * z + c;
    }
    complex double c_new = c - (z - target) / dc;
    double d2 = cabs2(c_new - c);
    if (cisfinite(c_new)) {
      c = c_new;
    } else {
      break;
    }
    if (d2 <= epsilon2) {
      break;
    }
  }
  ray->j = ray->j + 1;
  double d2 = cabs2(c - ray->c);
  if (d2 <= epsilon2) {
    ray->c = c;
    return m_converged;
  }
  if (cisfinite(c)) {
    ray->c = c;
    return m_stepped;
  } else {
    return m_failed;
  }
}

extern complex double m_d_exray_in_get(const m_d_exray_in *ray) {
  if (! ray) {
    return 0;
  }
  return ray->c;
}
