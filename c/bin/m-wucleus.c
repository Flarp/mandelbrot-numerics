#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision z-guess-re z-guess-im c-re c-im period maxsteps\n"
    , progname
    );
}

extern int main(int argc, char **argv) {
  if (argc != 8) {
    usage(argv[0]);
    return 1;
  }
  bool native = true;
  int bits = 0;
  if (! arg_precision(argv[1], &native, &bits)) { return 1; }
  if (native) {
    double zre = 0;
    double zim = 0;
    double cre = 0;
    double cim = 0;
    int period = 0;
    int maxsteps = 0;
    if (! arg_double(argv[2], &zre)) { return 1; }
    if (! arg_double(argv[3], &zim)) { return 1; }
    if (! arg_double(argv[4], &cre)) { return 1; }
    if (! arg_double(argv[5], &cim)) { return 1; }
    if (! arg_int(argv[6], &period)) { return 1; }
    if (! arg_int(argv[7], &maxsteps)) { return 1; }
    complex double z = 0;
    m_d_wucleus(&z, zre + I * zim, cre + I * cim, period, maxsteps);
    printf("%.16e %.16e\n", creal(z), cimag(z));
    return 0;
  } else {
    fprintf(stderr, "non-double precision not supported yet\n");
  }
  return 1;
}
