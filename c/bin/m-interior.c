#include <stdio.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

static void usage(const char *progname) {
  fprintf
    ( stderr
    , "usage: %s precision z-guess-re z-guess-im c-guess-re c-guess-im interior-r interior-t period maxsteps\n"
    , progname
    );
}

extern int main(int argc, char **argv) {
  if (argc != 10) {
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
    double ir = 0;
    double it = 0;
    int period = 0;
    int maxsteps = 0;
    if (! arg_double(argv[2], &zre)) { return 1; }
    if (! arg_double(argv[3], &zim)) { return 1; }
    if (! arg_double(argv[4], &cre)) { return 1; }
    if (! arg_double(argv[5], &cim)) { return 1; }
    if (! arg_double(argv[6], &ir)) { return 1; }
    if (! arg_double(argv[7], &it)) { return 1; }
    if (! arg_int(argv[8], &period)) { return 1; }
    if (! arg_int(argv[9], &maxsteps)) { return 1; }
    complex double z = 0;
    complex double c = 0;
    m_d_interior(&z, &c, zre + I * zim, cre + I * cim, ir * cexp(I * twopi * it), period, maxsteps);
    printf("%.16e %.16e %.16e %.16e\n", creal(z), cimag(z), creal(c), cimag(c));
    return 0;
  } else {
    mpc_t z, c, i, interior;
    mpfr_t p;
    mpc_init2(z, bits);
    mpc_init2(c, bits);
    mpc_init2(i, bits);
    mpc_init2(interior, bits);
    mpfr_init2(p, bits);
    int period = 0;
    int maxsteps = 0;
    if (! arg_mpc(argv[2], argv[3], z)) { return 1; }
    if (! arg_mpc(argv[4], argv[5], c)) { return 1; }
    if (! arg_mpc(argv[6], argv[7], i)) { return 1; }
    if (! arg_int(argv[8], &period)) { return 1; }
    if (! arg_int(argv[9], &maxsteps)) { return 1; }
    // interior = ir * cexp(I * twopi * it);
    mpfr_const_pi(p, MPFR_RNDN);
    mpfr_mul(mpc_imagref(i), mpc_imagref(i), p, MPFR_RNDN);
    mpfr_mul_2si(mpc_imagref(i), mpc_imagref(i), 1, MPFR_RNDN);
    mpfr_sin_cos(mpc_imagref(interior), mpc_realref(interior), mpc_imagref(i), MPFR_RNDN);
    mpc_mul_fr(interior, interior, mpc_realref(i), MPC_RNDNN);
    m_r_interior(z, c, z, c, interior, period, maxsteps);
    mpfr_printf("%Re %Re %Re %Re\n", mpc_realref(z), mpc_imagref(z), mpc_realref(c), mpc_imagref(c));
    mpc_clear(z);
    mpc_clear(c);
    mpc_clear(i);
    mpc_clear(interior);
    mpfr_clear(p);
    return 0;
  }
  return 1;
}
