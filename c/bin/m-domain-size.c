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
    double size = m_d_domain_size(nre + I * nim, period);
    printf("%.16e\n", size);
    return 0;
  } else {
    mpc_t n;
    int period = 0;
    mpc_init2(n, bits);
    if (! arg_mpc(argv[2], argv[3], n)) { return 1; }
    if (! arg_int(argv[4], &period)) { return 1; }
    mpfr_t size;
    mpfr_init2(size, bits);
    m_r_domain_size(size, n, period);
    mpfr_printf("%Re\n", size);
    mpfr_clear(size);
    mpc_clear(n);
    return 0;
  }
  return 1;
}
