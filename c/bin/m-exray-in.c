#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision angle sharpness maxsteps\n"
    , progname
    );
}

extern int main(int argc, char **argv) {
  if (argc != 5) {
    usage(argv[0]);
    return 1;
  }
  bool native = true;
  int bits = 0;
  if (! arg_precision(argv[1], &native, &bits)) { return 1; }
  mpq_t angle;
  mpq_init(angle);
  int sharpness = 0;
  int maxsteps = 0;
  if (! arg_rational(argv[2], angle)) { mpq_clear(angle); return 1; }
  if (! arg_int(argv[3], &sharpness)) { mpq_clear(angle); return 1; }
  if (! arg_int(argv[4], &maxsteps)) { mpq_clear(angle); return 1; }
  int retval = 0;
  if (native) {
    m_d_exray_in *ray = m_d_exray_in_new(angle, sharpness);
    if (! ray) { mpq_clear(angle); return 1; }
    for (int i = 0; i < maxsteps; ++i) {
      complex double c = m_d_exray_in_get(ray);
      printf("%.16e %.16e\n", creal(c), cimag(c));
      if (m_stepped != m_d_exray_in_step(ray)) {
        retval = 1;
        break;
      }
    }
    m_d_exray_in_delete(ray);
  } else {
    m_r_exray_in *ray = m_r_exray_in_new(angle, sharpness);
    if (! ray) { mpq_clear(angle); return 1; }
    mpc_t c;
    mpc_init2(c, 53);
    for (int i = 0; i < maxsteps; ++i) {
      m_r_exray_in_get(ray, c);
      mpfr_printf("%Re %Re\n", mpc_realref(c), mpc_imagref(c));
      if (m_stepped != m_r_exray_in_step(ray)) {
        retval = 1;
        break;
      }
      // prec
      mpfr_prec_t precr, preci, prec;
      mpc_get_prec2(&precr, &preci, c);
      prec = precr > preci ? precr : preci;
      if (prec > bits) {
        break;
      }
    }
    m_r_exray_in_delete(ray);
  }
  mpq_clear(angle);
  return retval;
}
