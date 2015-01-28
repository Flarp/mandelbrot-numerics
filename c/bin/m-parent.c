#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision nucleus-re nucleus-im period maxsteps\n"
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
    mpq_t angle;
    mpq_init(angle);
    complex double root = 0;
    complex double parent = 0;
    int p = m_d_parent(angle, &root, &parent, cre + I * cim, period, maxsteps);
    if (p > 0) {
      gmp_printf("%.18e %.18e %Qd %.18e %.18e %d\n", creal(root), cimag(root), angle, creal(parent), cimag(parent), p);
    } else if (p == 0) {
      printf("%.18e %.18e\n", creal(root), cimag(root));
    }
    mpq_clear(angle);
    return p < 0;
  } else {
    fprintf(stderr, "non-double precision not supported yet\n");
  }
  return 1;
}
