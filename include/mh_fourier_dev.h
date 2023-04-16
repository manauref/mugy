/* mugy: mh_fourier_dev.h
 *
 * Fourier type for the device.
 *
 */
#pragma once

#include "mh_userFLAGS.h"

#if (__NVCC__ && USE_GPU)

#include <mh_cuComplex_wrap.h>  /* For complex data types. */

#ifdef USE_SINGLE_PRECISION
typedef cuComplex mugy_cufourier_t;
#else
typedef cuDoubleComplex mugy_cufourier_t;
#endif

#endif
