#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision guess-re guess-im period maxsteps\n"
    , progname
    );
}

extern int main(int argc, char **argv) {
  if (argc != 6) {
    usage(argv[0]);
    return 1;
  }
  bool native = true;
  int bits = 0;
  if (! arg_precision(argv[1], &native, &bits)) { return 1; }
  if (native) {
    double cre = 0;
    double cim = 0;
    int period = 0;
    int maxsteps = 0;
    if (! arg_double(argv[2], &cre)) { return 1; }
    if (! arg_double(argv[3], &cim)) { return 1; }
    if (! arg_int(argv[4], &period)) { return 1; }
    if (! arg_int(argv[5], &maxsteps)) { return 1; }
    complex double c = 0;
    m_d_nucleus(&c, cre + I * cim, period, maxsteps);
    printf("%.16e %.16e\n", creal(c), cimag(c));
    return 0;
  } else {
    mpc_t c_guess;
    int period = 0;
    int maxsteps = 0;
    mpc_init2(c_guess, bits);
    if (! arg_mpc(argv[2], argv[3], c_guess)) { return 1; }
    if (! arg_int(argv[4], &period)) { return 1; }
    if (! arg_int(argv[5], &maxsteps)) { return 1; }
    mpc_t c_out;
    mpc_init2(c_out, bits);
    m_r_nucleus(c_out, c_guess, period, maxsteps);
    mpfr_printf("%Re %Re\n", mpc_realref(c_out), mpc_imagref(c_out));
    mpc_clear(c_out);
    mpc_clear(c_guess);
    return 0;
  }
  return 1;
}
