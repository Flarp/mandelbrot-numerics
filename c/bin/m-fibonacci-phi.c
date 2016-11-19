#include <stdio.h>
#include <gmp.h>
#include <mandelbrot-numerics.h>
#include "m-util.h"

int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;
  mpq_t a, b, c;
  mpq_init(a);
  mpq_init(b);
  mpq_init(c);
  mpq_set_si(a, 0, 1);
  mpq_set_si(b, 1, 1);
  double prev_size = 0.0/0.0;
  for (int i = 0; i < 24; ++i)
  {
    mpz_add(mpq_numref(c), mpq_numref(a), mpq_numref(b));
    mpz_add(mpq_denref(c), mpq_denref(a), mpq_denref(b));
    mpq_canonicalize(c);
    complex double z = 0, bond = 0, nucleus = 0;
    m_d_interior(&z, &bond, 0, 0, cexp(I * twopi * mpq_get_d(c)), 1, 64);
    int period = mpz_get_si(mpq_denref(c));
    m_d_nucleus(&nucleus, bond, period, 64);
    double size = cabs(m_d_size(nucleus, period));
    printf
      ( "%d/%d\t%.18e\t%.18e\n"
      , (int) mpz_get_si(mpq_numref(c))
      , period
      , size
      , sqrt(prev_size / size)
      );
    fflush(stdout);
    mpq_set(a, b);
    mpq_set(b, c);
    prev_size = size;
  }
  mpq_clear(a);
  mpq_clear(b);
  mpq_clear(c);
  return 0;
}
