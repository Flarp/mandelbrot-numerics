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
    m_shape shape = m_d_shape(nre + I * nim, period);
    switch (shape) {
      case m_cardioid: printf("cardioid\n"); return 0;
      case m_circle:   printf("circle\n");   return 0;
    }
  } else {
    fprintf(stderr, "non-double precision not supported yet\n");
  }
  return 1;
}
