#ifndef MANDELBROT_NUMERICS_H
#define MANDELBROT_NUMERICS_H 1

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

/* functions returning bool return true for success, false for failure */

struct m_d_mat2 {
  double _Complex a, b, c, d;
};
typedef struct m_d_mat2 m_d_mat2;

extern void m_d_mat2_set(m_d_mat2 *o, const m_d_mat2 *m);
extern void m_d_mat2_id(m_d_mat2 *o);
extern double _Complex m_d_mat2_tr(const m_d_mat2 *m);
extern double _Complex m_d_mat2_det(const m_d_mat2 *m);
extern void m_d_mat2_inv(m_d_mat2 *m1, const m_d_mat2 *m);
extern void m_d_mat2_mul(m_d_mat2 *o, const m_d_mat2 *l, const m_d_mat2 *r);
extern void m_d_mat2_diagonalize(m_d_mat2 *p, m_d_mat2 *d, m_d_mat2 *p1, const m_d_mat2 *m);
extern void m_d_mat2_moebius3(m_d_mat2 *m, double _Complex zero, double _Complex one, double _Complex infinity);

struct m_d_mat2_interp;
typedef struct m_d_mat2_interp m_d_mat2_interp;

extern m_d_mat2_interp *m_d_mat2_interp_new(void);
extern void m_d_mat2_interp_delete(m_d_mat2_interp *i);
extern void m_d_mat2_interp_init(m_d_mat2_interp *i, const m_d_mat2 *f, const m_d_mat2 *g);
extern void m_d_mat2_interp_do(m_d_mat2 *m, const m_d_mat2_interp *i, double t);

enum m_shape { m_cardioid, m_circle };
typedef enum m_shape m_shape;

enum m_newton { m_failed, m_stepped, m_converged };
typedef enum m_newton m_newton;

/* double precision: m_d_*()  */
/* functions taking non-const mpq_t use them for output */
/* functions taking pointers to double _Complex use them for output */

extern m_newton m_d_nucleus_step(double _Complex *c, double _Complex c_guess, int period);
extern m_newton m_d_nucleus(double _Complex *c, double _Complex c_guess, int period, int maxsteps);

extern m_newton m_d_misiurewicz_naive_step(double _Complex *c_out, double _Complex c_guess, int preperiod, int period);
extern m_newton m_d_misiurewicz_naive(double _Complex *c_out, double _Complex c_guess, int preperiod, int period, int mxsteps);
extern m_newton m_d_misiurewicz_step(double _Complex *c_out, double _Complex c_guess, int preperiod, int period);
extern m_newton m_d_misiurewicz(double _Complex *c_out, double _Complex c_guess, int preperiod, int period, int mxsteps);

extern m_newton m_d_wucleus_step(double _Complex *z, double _Complex z_guess, double _Complex c, int period);
extern m_newton m_d_wucleus(double _Complex *z, double _Complex z_guess, double _Complex c, int period, int maxsteps);

extern m_newton m_d_interior_step(double _Complex *z, double _Complex *c, double _Complex z_guess, double _Complex c_guess, double _Complex interior, int period);
extern m_newton m_d_interior(double _Complex *z, double _Complex *c, double _Complex z_guess, double _Complex c_guess, double _Complex interior, int period, int maxsteps);

extern int m_d_parent(mpq_t angle, double _Complex *root_out, double _Complex *parent_out, double _Complex nucleus, int period, int maxsteps);

extern bool m_d_interior_de(double *de_out, double _Complex *dz_out, double _Complex z, double _Complex c, int p, int steps);

extern double _Complex m_d_size(double _Complex nucleus, int period);
extern double m_d_domain_size(double _Complex nucleus, int period);
extern m_shape m_d_shape(double _Complex nucleus, int period);

struct m_d_exray_in;
typedef struct m_d_exray_in m_d_exray_in;
extern m_d_exray_in *m_d_exray_in_new(const mpq_t angle, int sharpness);
extern void m_d_exray_in_delete(m_d_exray_in *ray);
extern m_newton m_d_exray_in_step(m_d_exray_in *ray, int maxsteps);
extern double _Complex m_d_exray_in_get(const m_d_exray_in *ray);
extern double _Complex m_d_exray_in_do(const mpq_t angle, int sharpness, int maxsteps, int maxnewtonsteps);

struct m_d_exray_out;
typedef struct m_d_exray_out m_d_exray_out;
extern m_d_exray_out *m_d_exray_out_new(double _Complex c, int sharpness, int maxdwell);
extern void m_d_exray_out_delete(m_d_exray_out *ray);
extern m_newton m_d_exray_out_step(m_d_exray_out *ray);
extern bool m_d_exray_out_have_bit(const m_d_exray_out *ray);
extern bool m_d_exray_out_get_bit(const m_d_exray_out *ray);
extern double _Complex m_d_exray_out_get(const m_d_exray_out *ray);
extern char *m_d_exray_out_do(double _Complex c, int sharpness, int maxdwell);

struct m_d_box_period;
typedef struct m_d_box_period m_d_box_period;
extern m_d_box_period *m_d_box_period_new(double _Complex center, double radius);
extern void m_d_box_period_delete(m_d_box_period *box);
extern bool m_d_box_period_step(m_d_box_period *box);
extern bool m_d_box_period_have_period(const m_d_box_period *box);
extern int m_d_box_period_get_period(const m_d_box_period *box);
extern int m_d_box_period_do(double _Complex center, double radius, int maxperiod);

/* arbitrary precision: m_r_*() */
/* functions taking non-const mpq_t/mpfr_t/mpc_t use them for output */

extern m_newton m_r_nucleus_step_raw(mpc_t c_out, const mpc_t c_guess, int period, mpc_t z, mpc_t dc, mpc_t c_new, mpc_t d, mpfr_t d2, mpfr_t epsilon2);
extern m_newton m_r_nucleus_step(mpc_t c_out, const mpc_t c_guess, int period);
extern m_newton m_r_nucleus(mpc_t c_out, const mpc_t c_guess, int period, int maxsteps);

extern m_newton m_r_wucleus_step_raw(mpc_t z_out, const mpc_t z_guess, const mpc_t c, int period, mpc_t z, mpc_t dz, mpc_t z_new, mpc_t d, mpfr_t d2, mpfr_t epsilon2);
extern m_newton m_r_wucleus_step(mpc_t z_out, const mpc_t z_guess, const mpc_t c, int period);
extern m_newton m_r_wucleus(mpc_t z_out, const mpc_t z_guess, const mpc_t c, int period, int maxsteps);

extern m_newton m_r_interior_step_raw(mpc_t z_out, mpc_t c_out, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period, mpc_t z_new, mpc_t c_new, mpc_t c, mpc_t z, mpc_t dz, mpc_t dc, mpc_t dzdz, mpc_t dcdz, mpc_t det, mpc_t d, mpc_t dz1, mpfr_t d2z, mpfr_t d2c, mpfr_t epsilon2);
extern m_newton m_r_interior_step(mpc_t z_out, mpc_t c_out, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period);
extern m_newton m_r_interior(mpc_t z_out, mpc_t c_out, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period, int maxsteps);

extern int m_r_parent(mpq_t angle_out, mpc_t root_out, mpc_t parent_out, const mpc_t nucleus, int period, int maxsteps);

extern void m_r_size(mpc_t size, const mpc_t nucleus, int period);
extern int m_r_domain_size(mpfr_t size, const mpc_t nucleus, int period);
extern m_shape m_r_shape(const mpc_t nucleus, int period);

struct m_r_exray_in;
typedef struct m_r_exray_in m_r_exray_in;
extern m_r_exray_in *m_r_exray_in_new(const mpq_t angle, int sharpness);
extern void m_r_exray_in_delete(m_r_exray_in *ray);
extern m_newton m_r_exray_in_step(m_r_exray_in *ray, int maxsteps);
extern void m_r_exray_in_get(const m_r_exray_in *ray, mpc_t c);

/*
struct m_r_exray_out;
typedef struct m_r_exray_out m_r_exray_out;
extern m_r_exray_out *m_r_exray_out_new(const mpc_t c);
extern void m_r_exray_out_delete(m_r_exray_out *ray);
extern bool m_r_exray_out_step(m_r_exray_out *ray);
extern void m_r_exray_out_get(const m_r_exray_out *ray, mpc_t endpoint);
extern bool m_r_exray_out_have_bit(const m_r_exray_out *ray);
extern bool m_r_exray_out_get_bit(const m_r_exray_out *ray);
*/

struct m_r_box_period;
typedef struct m_r_box_period m_r_box_period;
extern m_r_box_period *m_r_box_period_new(const mpc_t center, const mpfr_t radius);
extern void m_r_box_period_delete(m_r_box_period *box);
extern bool m_r_box_period_step(m_r_box_period *box);
extern bool m_r_box_period_have_period(const m_r_box_period *box);
extern int m_r_box_period_get_period(const m_r_box_period *box);
extern int m_r_box_period_do(const mpc_t center, const mpfr_t radius, int maxperiod);

#endif
