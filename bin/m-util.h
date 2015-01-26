#ifndef M_UTIL_H
#define M_UTIL_H 1

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static inline bool arg_precision(const char *arg, bool *native, int *bits) {
  if (0 == strcmp("double", arg)) {
    *native = true;
    *bits = 53;
    return true;
  } else {
    char *check = 0;
    errno = 0;
    long int li = strtol(arg, &check, 10);
    bool valid = ! errno && arg != check && ! *check;
    if (valid && li > 0) {
      *native = false;
      *bits = li;
      return true;
    }
  }
  return false;
}

static inline bool arg_double(const char *arg, double *x) {
  char *check = 0;
  errno = 0;
  double d = strtod(arg, &check);
  if (! errno && arg != check && ! *check) {
    *x = d;
    return true;
  }
  return false;
}

static inline bool arg_int(const char *arg, int *x) {
  char *check = 0;
  errno = 0;
  long int li = strtol(arg, &check, 10);
  if (! errno && arg != check && ! *check) {
    *x = li;
    return true;
  }
  return false;
}

static inline bool arg_rational(const char *arg, mpq_t x) {
  int ok = mpq_set_str(x, arg, 10);
  mpq_canonicalize(x);
  return ok == 0;
}

static const double twopi = 6.283185307179586;

#endif
