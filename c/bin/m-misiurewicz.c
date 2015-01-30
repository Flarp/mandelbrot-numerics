#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision guess-re guess-im preperiod period maxsteps\n"
    , progname
    );
}

extern int main(int argc, char **argv) {
  if (argc != 7) {
    usage(argv[0]);
    return 1;
  }
  bool native = true;
  int bits = 0;
  if (! arg_precision(argv[1], &native, &bits)) { return 1; }
  if (native) {
    double cre = 0;
    double cim = 0;
    int preperiod = 0;
    int period = 0;
    int maxsteps = 0;
    if (! arg_double(argv[2], &cre)) { return 1; }
    if (! arg_double(argv[3], &cim)) { return 1; }
    if (! arg_int(argv[4], &preperiod)) { return 1; }
    if (! arg_int(argv[5], &period)) { return 1; }
    if (! arg_int(argv[6], &maxsteps)) { return 1; }
    complex double c = 0;
    m_d_misiurewicz(&c, cre + I * cim, preperiod, period, maxsteps);
    printf("%.16e %.16e\n", creal(c), cimag(c));
    return 0;
  } else {
    fprintf(stderr, "non-double precision not supported yet\n");
  }
  return 1;
}
