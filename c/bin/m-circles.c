// m-circles.c (C) 2016 Claude Heiland-Allen <claude@mathr.co.uk>
// demonstrate distortion of some circle-like components in the Mandelbrot set
// see: http://math.stackexchange.com/q/1857237/209286

#include <complex.h>
#include <stdio.h>

#include <mandelbrot-numerics.h>
// https://code.mathr.co.uk/mandelbrot-numerics
//
// m_d_interior(z_out, c_out, z_guess, c_guess, coordinate, period, steps);
// Find a point with given interior 'coordinate' in a hyperbolic component of
// 'period' (magnitude 1 is boundary).
// Uses Newton's method with at most 'steps'.
//
// m_d_nucleus(c_out, c_guess, period, steps);
// Find the nucleus of a hyperbolic component with given 'period'.
// Uses Newton's method with at most 'steps'.
//
// m_d_size(nucleus, period);
// Find a size estimate for a hyperbolic component.
// For circular-like components, approximates the diameter.

#define twopi 6.283185307179586

int main(int argc, char **argv) {
  (void) argc;
  (void) argv;

  // initial component
  int period = 2;
  double size = 0.5;
  complex double nucleus = -1;
  complex double leftpoint = -1.25;
  complex double rightpoint = -0.75;

  // consider components heading towards -inf on the period-doubling cascade
  for (int depth = 0; depth < 12; ++depth) {

    // geometric midpoint might differ from nucleus in general
    complex double midpoint = (leftpoint + rightpoint) / 2;
    double radius = cabs(rightpoint - midpoint);

    // trace boundaries of components
    complex double z = nucleus, c = nucleus;
    for (int t = 0; t < 360; ++t) {

      // m_d_interior is unstable at root point (theta = 0), hence adding 0.5deg
      double theta = twopi * (t + 0.5) / 360;
      m_d_interior(&z, &c, z, c, cexp(I * theta), period, 64);

      // compute the geometric distance to the midpoint
      // divide by the geometric radius to normalize output range
      double distance = cabs(c - midpoint) / radius;

      // output with full precision
      printf("%.18g\n", distance);
    }

    // separator between curves for gnuplot
    printf("\n\n");

    // advance to the next component
    period *= 2;

    // the next component is approximately 1/4 the size
    complex double estimate = nucleus - 1.25 * (size / 2);
    m_d_nucleus(&nucleus, estimate, period, 64);
    size = cabs(m_d_size(nucleus, period));

    // update edges - can't find rightpoint with m_d_interior due to instability
    rightpoint = leftpoint;
    m_d_interior(&z, &leftpoint, nucleus, nucleus, -1, period, 64);

  }

  return 0;
}
