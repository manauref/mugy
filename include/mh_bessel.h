/* mugy: mh_bessel.h
 *
 * Functions that compute the Bessel functions needed.
 *
 */
#include "mh_macros.h"
#include <gsl/gsl_sf_bessel.h>

#ifdef USE_SINGLE_PRECISION

static inline real mugy_bessel_I0_scaled(real x) {
  return (real) gsl_sf_bessel_I0_scaled((double) x);
}

static inline real mugy_bessel_I1_scaled(real x) {
  return (real) gsl_sf_bessel_I1_scaled((double) x);
}

#else
// ............ DOUBLE PRECISION ............ //

static inline real mugy_bessel_I0_scaled(real x) {
  return gsl_sf_bessel_I0_scaled(x);
}

static inline real mugy_bessel_I1_scaled(real x) {
  return gsl_sf_bessel_I1_scaled(x);
}

#endif
