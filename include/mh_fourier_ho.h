/* mugy: mh_fourier_ho.h
 *
 * Fourier type for the host.
 *
 */
#pragma once

#include "mh_userFLAGS.h"
#include <complex.h>  /* For complex data types. */

#ifdef USE_SINGLE_PRECISION
typedef float complex fourier;
#else
typedef double complex fourier;
#endif
