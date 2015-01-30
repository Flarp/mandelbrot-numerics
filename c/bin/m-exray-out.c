#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision c-re c-im sharpness maxdwell preperiod period\n"
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
    double cre = 0;
    double cim = 0;
    int sharpness = 0;
    int maxdwell = 0;
    int preperiod = 0;
    int period = 0;
    if (! arg_double(argv[2], &cre)) { return 1; }
    if (! arg_double(argv[3], &cim)) { return 1; }
    if (! arg_int(argv[4], &sharpness)) { return 1; }
    if (! arg_int(argv[5], &maxdwell)) { return 1; }
    if (! arg_int(argv[6], &preperiod)) { return 1; }
    if (! arg_int(argv[7], &period)) { return 1; }
    if (! (preperiod >= 0)) { return 1; }
    if (! (period >= 0)) { return 1; }
    if (preperiod > 0 && ! (period > 0)) { return 1; }
    m_d_exray_out *ray = m_d_exray_out_new(cre + I * cim, sharpness, maxdwell);
    if (! ray) { return 1; }
    char *bits = malloc(maxdwell + 2);
    if (! bits) { m_d_exray_out_delete(ray); return 1; }
    int n = 0;
    do {
      complex double c = m_d_exray_out_get(ray);
      printf("%.16e %.16e", creal(c), cimag(c));
      if (m_d_exray_out_have_bit(ray)) {
        int bit = m_d_exray_out_get_bit(ray);
        bits[n++] = '0' + bit;
        printf(" %d\n", bit);
      } else {
        printf("\n");
      }
    } while (m_stepped == m_d_exray_out_step(ray));
    bits[n] = 0;
    m_d_exray_out_delete(ray);
    if (preperiod + period > 0) {
      printf("\n.");
      for (int i = 0; i < preperiod && 0 <= n - 1 - i; ++i) {
        putchar(bits[n - 1 - i]);
      }
      putchar('(');
      for (int i = preperiod; i < preperiod + period && 0 <= n - 1 - i; ++i) {
        putchar(bits[n - 1 - i]);
      }
      putchar(')');
      putchar('\n');
    }
    free(bits);
    return 0;
  } else {
    fprintf(stderr, "non-double precision not supported yet\n");
  }
  return 1;
}
