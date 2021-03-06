#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision center-re center-im radius maxperiod\n"
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
    double radius = 0;
    int maxperiod = 0;
    if (! arg_double(argv[2], &cre)) { return 1; }
    if (! arg_double(argv[3], &cim)) { return 1; }
    if (! arg_double(argv[4], &radius)) { return 1; }
    if (! arg_int(argv[5], &maxperiod)) { return 1; }
    int period = m_d_box_period_do(cre + I * cim, radius, maxperiod);
    if (period > 0) {
      printf("%d\n", period);
      return 0;
    }
  } else {
    mpc_t center;
    mpfr_t radius;
    int maxperiod = 0;
    mpc_init2(center, bits);
    mpfr_init2(radius, bits);
    if (! arg_mpc(argv[2], argv[3], center)) { return 1; }
    if (! arg_mpfr(argv[4], radius)) { return 1; }
    if (! arg_int(argv[5], &maxperiod)) { return 1; }
    int period = m_r_box_period_do(center, radius, maxperiod);
    if (period > 0) {
      printf("%d\n", period);
    }
    mpc_clear(center);
    mpfr_clear(radius);
    return period <= 0;
  }
  return 1;
}
