#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision nucleus-re nucleus-im period\n"
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
  if (native) {
    double nre = 0;
    double nim = 0;
    int period = 0;
    if (! arg_double(argv[2], &nre)) { return 1; }
    if (! arg_double(argv[3], &nim)) { return 1; }
    if (! arg_int(argv[4], &period)) { return 1; }
    double _Complex s = m_d_shape_estimate(nre + I * nim, period);
    m_shape shape = m_d_shape_discriminant(s);
    switch (shape) {
      case m_cardioid: printf("cardioid %.16e %.16e\n", creal(s), cimag(s)); return 0;
      case m_circle:   printf("circle %.16e %.16e\n",   creal(s), cimag(s)); return 0;
    }
  } else {
    mpc_t n;
    int period = 0;
    mpc_init2(n, bits);
    if (! arg_mpc(argv[2], argv[3], n)) { return 1; }
    if (! arg_int(argv[4], &period)) { return 1; }
    mpc_t s;
    mpc_init2(s, 53);
    m_r_shape_estimate(s, n, period);
    m_shape shape = m_r_shape_discriminant(s);
    mpc_clear(n);
    switch (shape) {
      case m_cardioid: mpfr_printf("cardioid %Re %Re\n", mpc_realref(s), mpc_imagref(s)); break;
      case m_circle:   mpfr_printf("circle %Re %Re\n"  , mpc_realref(s), mpc_imagref(s)); break;
    }
    mpc_clear(s);
  }
  return 1;
}
