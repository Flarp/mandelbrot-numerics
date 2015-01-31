#include <mandelbrot-numerics.h>
#include "m_d_util.h"

struct m_r_exray_in {
  mpq_t angle;
  mpq_t one;
  int sharpness;
  double er;
  mpfr_prec_t prec;
  mpc_t c0, c, c_new, target, z, dc, d;
  mpfr_t d2, epsilon2, epsilon256;
  int j;
  int k;
};

extern m_r_exray_in *m_r_exray_in_new(const mpq_t angle, int sharpness) {
  m_r_exray_in *ray = malloc(sizeof(*ray));
  mpq_init(ray->angle);
  mpq_set(ray->angle, angle);
  mpq_init(ray->one);
  mpq_set_ui(ray->one, 1, 1);
  ray->sharpness = sharpness;
  ray->er = 65536.0;
  ray->prec = 53;
  mpc_init2(ray->c0, ray->prec);
  mpc_init2(ray->c, ray->prec);
  mpc_init2(ray->c_new, ray->prec);
  mpc_init2(ray->target, ray->prec);
  mpc_init2(ray->z, ray->prec);
  mpc_init2(ray->dc, ray->prec);
  mpc_init2(ray->d, ray->prec);
  mpfr_init2(ray->d2, ray->prec);
  mpfr_init2(ray->epsilon2, ray->prec);
  mpfr_init2(ray->epsilon256, ray->prec);
  double a = twopi * mpq_get_d(ray->angle);
  // c0 = er * (cos(a) + I * sin(a));
  mpc_set_d_d(ray->c0, ray->er * cos(a), ray->er * sin(a), MPC_RNDNN);
  ray->k = 0;
  ray->j = 0;
  // epsilon
  mpfr_set_si(ray->epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(ray->epsilon2);
  mpfr_sub_si(ray->epsilon2, ray->epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(ray->epsilon2, ray->epsilon2, MPFR_RNDN);
  mpfr_set_si(ray->epsilon256, 256, MPFR_RNDN);
  mpfr_nextabove(ray->epsilon256);
  mpfr_sub_si(ray->epsilon256, ray->epsilon256, 256, MPFR_RNDN);
  mpfr_sqr(ray->epsilon256, ray->epsilon256, MPFR_RNDN);
  return ray;
}

extern void m_r_exray_in_delete(m_r_exray_in *ray) {
  if (! ray) {
    return;
  }
  mpq_clear(ray->angle);
  mpq_clear(ray->one);
  mpc_clear(ray->c0);
  mpc_clear(ray->c);
  mpc_clear(ray->c_new);
  mpc_clear(ray->target);
  mpc_clear(ray->z);
  mpc_clear(ray->dc);
  mpc_clear(ray->d);
  mpfr_clear(ray->d2);
  mpfr_clear(ray->epsilon2);
  mpfr_clear(ray->epsilon256);
  free(ray);
}

static void m_r_exray_in_bump_prec(m_r_exray_in *ray) {
  ray->prec += 16;
  mpfr_prec_round(mpc_realref(ray->c0), ray->prec, MPFR_RNDN);
  mpfr_prec_round(mpc_imagref(ray->c0), ray->prec, MPFR_RNDN);
  mpc_set_prec(ray->c, ray->prec);
  mpc_set_prec(ray->c_new, ray->prec);
  mpc_set_prec(ray->target, ray->prec);
  mpc_set_prec(ray->z, ray->prec);
  mpc_set_prec(ray->dc, ray->prec);
  mpc_set_prec(ray->d, ray->prec);
  mpfr_set_prec(ray->d2, ray->prec);
  mpfr_set_prec(ray->epsilon2, ray->prec);
  mpfr_set_prec(ray->epsilon256, ray->prec);
  // epsilon
  mpfr_set_si(ray->epsilon2, 2, MPFR_RNDN);
  mpfr_nextabove(ray->epsilon2);
  mpfr_sub_si(ray->epsilon2, ray->epsilon2, 2, MPFR_RNDN);
  mpfr_sqr(ray->epsilon2, ray->epsilon2, MPFR_RNDN);
  mpfr_set_si(ray->epsilon256, 256, MPFR_RNDN);
  mpfr_nextabove(ray->epsilon256);
  mpfr_sub_si(ray->epsilon256, ray->epsilon256, 256, MPFR_RNDN);
  mpfr_sqr(ray->epsilon256, ray->epsilon256, MPFR_RNDN);
}

extern m_newton m_r_exray_in_step(m_r_exray_in *ray) {
  if (ray->j >= ray->sharpness) {
    mpq_mul_2exp(ray->angle, ray->angle, 1);
    if (mpq_cmp_ui(ray->angle, 1, 1) >= 0) {
      mpq_sub(ray->angle, ray->angle, ray->one);
    }
    ray->k = ray->k + 1;
    ray->j = 0;
  }
  int bumps = 0;
restart:
  if (! (bumps < 64)) {
    return m_failed;
  }
  double r = pow(ray->er, pow(0.5, (ray->j + 0.5) / ray->sharpness));
  double a = twopi * mpq_get_d(ray->angle);
  mpc_set_d_d(ray->target, r * cos(a), r * sin(a), MPC_RNDNN);
  // c = c0;
  mpc_set(ray->c, ray->c0, MPC_RNDNN);
  for (int i = 0; i < 64; ++i) { // FIXME arbitrary limit
    // z = 0; dc = 0;
    mpc_set_si(ray->z, 0, MPC_RNDNN);
    mpc_set_si(ray->dc, 0, MPC_RNDNN);
    for (int p = 0; p <= ray->k; ++p) {
      // dc = 2 * z * dc + 1;
      mpc_mul(ray->dc, ray->z, ray->dc, MPC_RNDNN);
      mpc_mul_2si(ray->dc, ray->dc, 1, MPC_RNDNN);
      mpc_add_ui(ray->dc, ray->dc, 1, MPC_RNDNN);
      // z = z * z + c;
      mpc_sqr(ray->z, ray->z, MPC_RNDNN);
      mpc_add(ray->z, ray->z, ray->c, MPC_RNDNN);
    }
    // c_new = c - (z - target) / dc;
    mpc_sub(ray->z, ray->z, ray->target, MPC_RNDNN);
    mpc_div(ray->z, ray->z, ray->dc, MPC_RNDNN);
    mpc_sub(ray->c_new, ray->c, ray->z, MPC_RNDNN);
    if (mpfr_regular_p(mpc_realref(ray->c_new)) && mpfr_regular_p(mpc_imagref(ray->c_new))) {
      // d2 = norm(c_new - c);
      mpc_sub(ray->d, ray->c_new, ray->c, MPC_RNDNN);
      mpc_norm(ray->d2, ray->d, MPFR_RNDN);
      // c = c_new;
      mpc_set(ray->c, ray->c_new, MPC_RNDNN);
      if (mpfr_lessequal_p(ray->d2, ray->epsilon2)) {
        break;
      }
    } else {
      break;
    }
  }
  ray->j = ray->j + 1;
  if (mpfr_regular_p(mpc_realref(ray->c)) && mpfr_regular_p(mpc_imagref(ray->c))) {
    // d2 = norm(c - c0);
    mpc_sub(ray->d, ray->c, ray->c0, MPC_RNDNN);
    mpc_norm(ray->d2, ray->d, MPFR_RNDN);
    if (mpfr_lessequal_p(ray->d2, ray->epsilon256)) {
      m_r_exray_in_bump_prec(ray);
      bumps++;
      goto restart;
    } else {
      // c0 = c;
      mpc_set(ray->c0, ray->c, MPC_RNDNN);
      return m_stepped;
    }
  } else {
    return m_failed;
  }
}

extern void m_r_exray_in_get(const m_r_exray_in *ray, mpc_t c) {
  if (! ray) {
    return;
  }
  mpc_set_prec(c, ray->prec);
  mpc_set(c, ray->c0, MPC_RNDNN);
}
