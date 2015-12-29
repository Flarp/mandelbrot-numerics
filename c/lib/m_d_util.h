#ifndef M_D_UTIL_H
#define M_D_UTIL_H 1

static inline int sgn(double z) {
  if (z > 0) { return  1; }
  if (z < 0) { return -1; }
  return 0;
}

static inline bool odd(int a) {
  return a & 1;
}

static inline double cabs2(double _Complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

static inline bool cisfinite(double _Complex z) {
  return isfinite(creal(z)) && isfinite(cimag(z));
}

static const double pi = 3.141592653589793;
static const double twopi = 6.283185307179586;

// last . takeWhile (\x -> 2 /= 2 + x) . iterate (/2) $ 1 :: Double
static const double epsilon = 4.440892098500626e-16;

// epsilon^2
static const double epsilon2 = 1.9721522630525295e-31;

#endif
