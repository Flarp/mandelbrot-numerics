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

enum m_shape { m_cardioid, m_circle };
typedef enum m_shape m_shape;

enum m_newton { m_failed, m_stepped, m_converged };
typedef enum m_newton m_newton;

/* double precision: m_d_*()  */
/* functions taking non-const mpq_t use them for output */
/* functions taking pointers to complex double use them for output */

extern m_newton m_d_nucleus(complex double *c, complex double c_guess, int period);
extern m_newton m_d_wucleus(complex double *z, complex double z_guess, complex double c, int period);
extern m_newton m_d_interior(complex double *z, complex double *c, complex double z_guess, complex double c_guess, complex double interior, int period);
extern int m_d_parent(mpq_t angle, complex double *root_out, complex double *parent_out, complex double nucleus, int period, int maxsteps);
extern complex double m_d_size(complex double nucleus, int period);
extern double m_d_domain_size(complex double nucleus, int period);
extern m_shape m_d_shape(complex double nucleus, int period);

struct m_d_exray_in;
typedef struct m_d_exray_in m_d_exray_in;
extern m_d_exray_in *m_d_exray_in_new(const mpq_t angle, int sharpness);
extern void m_d_exray_in_delete(m_d_exray_in *ray);
extern m_newton m_d_exray_in_step(m_d_exray_in *ray);
extern complex double m_d_exray_in_get(const m_d_exray_in *ray);

struct m_d_exray_out;
typedef struct m_d_exray_out m_d_exray_out;
extern m_d_exray_out *m_d_exray_out_new(complex double c, int sharpness, int maxdwell);
extern void m_d_exray_out_delete(m_d_exray_out *ray);
extern m_newton m_d_exray_out_step(m_d_exray_out *ray);
extern bool m_d_exray_out_have_bit(const m_d_exray_out *ray);
extern bool m_d_exray_out_get_bit(const m_d_exray_out *ray);
extern complex double m_d_exray_out_get(const m_d_exray_out *ray);

struct m_d_box_period;
typedef struct m_d_box_period m_d_box_period;
extern m_d_box_period *m_d_box_period_new(complex double center, double radius);
extern void m_d_box_period_delete(m_d_box_period *box);
extern bool m_d_box_period_step(m_d_box_period *box);
extern bool m_d_box_period_have_period(const m_d_box_period *box);
extern int m_d_box_period_get_period(const m_d_box_period *box);

/* arbitrary precision: m_r_*() */
/* functions taking non-const mpq_t/mpfr_t/mpc_t use them for output */

extern bool m_r_nucleus(mpc_t c, const mpc_t c_guess, int period, int maxsteps);
extern bool m_r_wucleus(mpc_t z, mpc_t dz, const mpc_t z_guess, const mpc_t c, int period, int maxsteps);
extern bool m_r_interior(mpc_t z, mpc_t c, const mpc_t z_guess, const mpc_t c_guess, const mpc_t interior, int period, int maxsteps);
extern bool m_r_parent(mpq_t angle, mpc_t root, mpc_t parent, const mpc_t nucleus, int period, int maxsteps);
extern void m_r_size(mpc_t size, const mpc_t nucleus, int period);
extern void m_r_domain_size(mpfr_t size, const mpc_t nucleus, int period);
extern m_shape m_r_shape(const mpc_t nucleus, int period);

struct m_r_exray_in;
typedef struct m_r_exray_in m_r_exray_in;
extern m_r_exray_in *m_r_exray_in_new(const mpq_t angle);
extern void m_r_exray_in_delete(m_r_exray_in *ray);
extern bool m_r_exray_in_step(m_r_exray_in *ray);
extern void m_r_exray_in_get(const m_r_exray_in *ray, mpc_t endpoint);

struct m_r_exray_out;
typedef struct m_r_exray_out m_r_exray_out;
extern m_r_exray_out *m_r_exray_out_new(const mpc_t c);
extern void m_r_exray_out_delete(m_r_exray_out *ray);
extern bool m_r_exray_out_step(m_r_exray_out *ray);
extern void m_r_exray_out_get(const m_r_exray_out *ray, mpc_t endpoint);
extern bool m_r_exray_out_have_bit(const m_r_exray_out *ray);
extern bool m_r_exray_out_get_bit(const m_r_exray_out *ray);

struct m_r_box_period;
typedef struct m_r_box_period m_r_box_period;
extern m_r_box_period *m_r_box_period_new(const mpc_t center, const mpfr_t radius);
extern void m_r_box_period_delete(m_r_box_period *box);
extern bool m_r_box_period_step(m_r_box_period *box);
extern bool m_r_box_period_have_period(const m_r_box_period *box);
extern int m_r_box_period_get_period(const m_r_box_period *box);

#endif
