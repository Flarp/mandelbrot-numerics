#include <mandelbrot-numerics.h>
#include "m_d_util.h"

static void cross(mpfr_t out, const mpc_t a, const mpc_t b, mpfr_t t) {
  mpfr_mul(out, mpc_imagref(a), mpc_realref(b), MPFR_RNDN);
  mpfr_mul(t, mpc_realref(a), mpc_imagref(b), MPFR_RNDN);
  mpfr_sub(out, out, t, MPFR_RNDN);
}

static bool crosses_positive_real_axis(const mpc_t a, const mpc_t b, mpc_t d, mpfr_t t1, mpfr_t t2) {
  if (mpfr_sgn(mpc_imagref(a)) != mpfr_sgn(mpc_imagref(b))) {
    // d = b - a;
    mpc_sub(d, b, a, MPC_RNDNN);
    // s = sgn(cimag(d));
    int s = mpfr_sgn(mpc_imagref(d));
    // t = sgn(cross(d, a));
    cross(t1, d, a, t2);
    int t = mpfr_sgn(t1);
    return s == t;
  }
  return false;
}

static bool surrounds_origin(const mpc_t a, const mpc_t b, const mpc_t c, const mpc_t d, mpc_t e, mpfr_t t1, mpfr_t t2) {
  return odd
    ( crosses_positive_real_axis(a, b, e, t1, t2)
    + crosses_positive_real_axis(b, c, e, t1, t2)
    + crosses_positive_real_axis(c, d, e, t1, t2)
    + crosses_positive_real_axis(d, a, e, t1, t2)
    );
}

struct m_r_box_period {
  mpc_t c[4];
  mpc_t z[4];
  int p;
  mpc_t e;
  mpfr_t t1;
  mpfr_t t2;
};

extern m_r_box_period *m_r_box_period_new(const mpc_t center, const mpfr_t radius) {
  m_r_box_period *box = malloc(sizeof(*box));
  if (! box) {
    return 0;
  }
  // prec
  mpfr_prec_t precr, preci, prec;
  mpc_get_prec2(&precr, &preci, center);
  prec = precr > preci ? precr : preci;
  // init
  for (int i = 0; i < 4; ++i) {
    mpc_init2(box->c[i], prec);
    mpc_init2(box->z[i], prec);
  }
  mpc_init2(box->e, prec);
  mpfr_init2(box->t1, prec);
  mpfr_init2(box->t2, prec);
  // box
  mpfr_sub(mpc_realref(box->c[0]), mpc_realref(center), radius, MPFR_RNDN);
  mpfr_sub(mpc_imagref(box->c[0]), mpc_imagref(center), radius, MPFR_RNDN);
  mpfr_add(mpc_realref(box->c[1]), mpc_realref(center), radius, MPFR_RNDN);
  mpfr_sub(mpc_imagref(box->c[1]), mpc_imagref(center), radius, MPFR_RNDN);
  mpfr_add(mpc_realref(box->c[2]), mpc_realref(center), radius, MPFR_RNDN);
  mpfr_add(mpc_imagref(box->c[2]), mpc_imagref(center), radius, MPFR_RNDN);
  mpfr_sub(mpc_realref(box->c[3]), mpc_realref(center), radius, MPFR_RNDN);
  mpfr_add(mpc_imagref(box->c[3]), mpc_imagref(center), radius, MPFR_RNDN);
  mpc_set(box->z[0], box->c[0], MPC_RNDNN);
  mpc_set(box->z[1], box->c[1], MPC_RNDNN);
  mpc_set(box->z[2], box->c[2], MPC_RNDNN);
  mpc_set(box->z[3], box->c[3], MPC_RNDNN);
  box->p = 1;
  return box;
}

extern void m_r_box_period_delete(m_r_box_period *box) {
  if (box) {
    for (int i = 0; i < 4; ++i) {
      mpc_clear(box->c[i]);
      mpc_clear(box->z[i]);
    }
    mpc_clear(box->e);
    mpfr_clear(box->t1);
    mpfr_clear(box->t2);
    free(box);
  }
}

extern bool m_r_box_period_step(m_r_box_period *box) {
  if (! box) {
    return false;
  }
  bool ok = true;
  for (int i = 0; i < 4; ++i) {
    // box->z[i] = box->z[i] * box->z[i] + box->c[i];
    mpc_sqr(box->z[i], box->z[i], MPC_RNDNN);
    mpc_add(box->z[i], box->z[i], box->c[i], MPC_RNDNN);
    ok = ok && mpfr_number_p(mpc_realref(box->z[i])) && mpfr_number_p(mpc_imagref(box->z[i]));
  }
  box->p = box->p + 1;
  return ok;
}

extern bool m_r_box_period_have_period(const m_r_box_period *box) {
  if (! box) {
    return true;
  }
  m_r_box_period *ubox = (m_r_box_period *) box; // const cast
  return surrounds_origin(box->z[0], box->z[1], box->z[2], box->z[3], ubox->e, ubox->t1, ubox->t2);
}

extern int m_r_box_period_get_period(const m_r_box_period *box) {
  if (! box) {
    return 0;
  }
  return box->p;
}

extern int m_r_box_period_do(const mpc_t center, const mpfr_t radius, int maxperiod) {
  m_r_box_period *box = m_r_box_period_new(center, radius);
  if (! box) {
    return 0;
  }
  int period = 0;
  for (int i = 0; i < maxperiod; ++i) {
    if (m_r_box_period_have_period(box)) {
      period = m_r_box_period_get_period(box);
      break;
    }
    if (! m_r_box_period_step(box)) {
      break;
    }
  }
  m_r_box_period_delete(box);
  return period;
}
