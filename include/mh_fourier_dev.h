/* mugy: mh_fourier_dev.h
 *
 * Fourier type for the device.
 *
 */
#ifndef MUGY_FOURIER_DEV
#define MUGY_FOURIER_DEV

#include "mh_userFLAGS.h"

#if (__NVCC__ && USE_GPU)

#include <cuComplex.h>  /* For complex data types. */

#if USE_SINGLE_PRECISION
typedef cuComplex cufourier;
#else
typedef cuDoubleComplex cufourier;
#endif

#endif

#endif
