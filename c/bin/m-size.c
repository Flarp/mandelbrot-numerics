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
    complex double size = m_d_size(nre + I * nim, period);
    printf("%.16e %.16e\n", cabs(size), carg(size));
    return 0;
  } else {
    mpc_t n;
    int period = 0;
    mpc_init2(n, bits);
    if (! arg_mpc(argv[2], argv[3], n)) { return 1; }
    if (! arg_int(argv[4], &period)) { return 1; }
    mpc_t size;
    mpc_init2(size, 53);
    m_r_size(size, n, period);
    mpfr_t r, t;
    mpfr_init2(r, 53);
    mpfr_init2(t, 53);
    mpc_abs(r, size, MPFR_RNDN);
    mpc_arg(t, size, MPFR_RNDN);
    mpfr_printf("%Re %Re\n", r, t);
    mpfr_clear(r);
    mpfr_clear(t);
    mpc_clear(size);
    mpc_clear(n);
    return 0;
  }
  return 1;
}
