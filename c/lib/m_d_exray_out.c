#include <mandelbrot-numerics.h>
#include "m_d_util.h"

static double dwell(double loger2, int n, double zmag2) {
  return n - log2(log(zmag2) / loger2);
}

struct m_d_exray_out {
  int sharpness;
  double er;
  double er2;
  double loger2;
  double _Complex c;
  double _Complex z;
  double d;
  int n;
  int bit;
};

extern m_d_exray_out *m_d_exray_out_new(double _Complex c, int sharpness, int maxdwell) {
  double er = 65536;
  double er2 = er * er;
  int n = 0;
  double _Complex z = 0;
  for (int i = 0; i < maxdwell; ++i) {
    n = n + 1;
    z = z * z + c;
    if (cabs2(z) > er2) {
      break;
    }
  }
  if (! (cabs2(z) > er2)) {
    return 0;
  }
  m_d_exray_out *ray = malloc(sizeof(*ray));
  if (! ray) {
    return 0;
  }
  ray->sharpness = sharpness;
  ray->er = er;
  ray->er2 = er2;
  ray->loger2 = log(er2);
  ray->c = c;
  ray->z = z;
  ray->d = dwell(ray->loger2, n, cabs2(z));
  ray->n = n;
  ray->bit = -1;
  return ray;
}

extern void m_d_exray_out_delete(m_d_exray_out *ray) {
  if (ray) {
    free(ray);
  }
}

extern m_newton m_d_exray_out_step(m_d_exray_out *ray) {
  if (! ray) {
    return m_failed;
  }
  ray->bit = -1;
  ray->d -= 1.0 / ray->sharpness;
  if (ray->d <= 0) {
    return m_converged;
  }
  int m = ceil(ray->d);
  double r = pow(ray->er, pow(2, m - ray->d));
  double a = carg(ray->z) / twopi;
  double t = a - floor(a);
  if (m == ray->n) {
    double _Complex k = r * cexp(I * twopi *  t);
    double _Complex c = ray->c;
    double _Complex z = 0;
    for (int i = 0; i < 64; ++i) { // FIXME arbitrary limit
      double _Complex dc = 0;
      z = 0;
      for (int p = 0; p < m; ++p) {
        dc = 2 * z * dc + 1;
        z = z * z + c;
      }
      double _Complex c_new = c - (z - k) / dc;
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
    if (cisfinite(c)) {
      ray->c = c;
      ray->z = z;
      ray->d = dwell(ray->loger2, m, cabs2(z));
      return m_stepped;
    }
    return m_failed;
  } else {
    double _Complex k[2] = { r * cexp(I * pi * t), r * cexp(I * pi * (t + 1)) };
    double _Complex c[2] = { ray->c, ray->c };
    double _Complex z[2] = { 0, 0 };
    double d2[2];
    double e2[2];
    for (int i = 0; i < 64; ++i) { // FIXME arbitrary limit
      z[0] = 0;
      z[1] = 0;
      double _Complex dc[2] = { 0, 0 };
      for (int p = 0; p < m; ++p) {
        for (int w = 0; w < 2; ++w) {
          dc[w] = 2 * z[w] * dc[w] + 1;
          z[w] = z[w] * z[w] + c[w];
        }
      }
      double _Complex c_new[2];
      for (int w = 0; w < 2; ++w) {
        c_new[w] = c[w] - (z[w] - k[w]) / dc[w];
        e2[w] = cabs2(c_new[w] - c[w]);
        d2[w] = cabs2(c_new[w] - ray->c);
        c[w] = c_new[w];
      }
      if (! (e2[0] > epsilon2 && e2[1] > epsilon2)) {
        break;
      }
    }
    if ((cisfinite(c[0]) && cisfinite(c[1]) && d2[0] <= d2[1])
     || (cisfinite(c[0]) && ! cisfinite(c[1]))) {
      ray->bit = 0;
      ray->c = c[0];
      ray->z = z[0];
      ray->n = m;
      ray->d = dwell(ray->loger2, m, cabs2(z[0]));
      return m_stepped;
    }
    if ((cisfinite(c[0]) && cisfinite(c[1]) && d2[1] <= d2[0])
     || (! cisfinite(c[0]) && cisfinite(c[1]))) {
      ray->bit = 1;
      ray->c = c[1];
      ray->z = z[1];
      ray->n = m;
      ray->d = dwell(ray->loger2, m, cabs2(z[1]));
      return m_stepped;
    }
    return m_failed;
  }
}

extern bool m_d_exray_out_have_bit(const m_d_exray_out * ray) {
  if (! ray) {
    return false;
  }
  return 0 <= ray->bit;
}

extern bool m_d_exray_out_get_bit(const m_d_exray_out *ray) {
  if (! ray) {
    return false;
  }
  return ray->bit;
}

extern double _Complex m_d_exray_out_get(const m_d_exray_out *ray) {
  if (! ray) {
    return 0;
  }
  return ray->c;
}

extern char *m_d_exray_out_do(double _Complex c, int sharpness, int maxdwell) {
  m_d_exray_out *ray = m_d_exray_out_new(c, sharpness, maxdwell);
  if (! ray) {
    return 0;
  }
  char *bits = malloc(maxdwell + 2);
  if (! bits) {
    m_d_exray_out_delete(ray);
    return 0;
  }
  int n = 0;
  while (n <= maxdwell) {
    if (m_d_exray_out_have_bit(ray)) {
      bits[n++] = '0' + m_d_exray_out_get_bit(ray);
    }
    if (m_stepped != m_d_exray_out_step(ray)) {
      break;
    }
  }
  bits[n] = 0;
  m_d_exray_out_delete(ray);
  return bits;
}
