/* mugy: mh_data_dev.h

   Data types and some macros used in mugy device (GPU) operations.
*/
#ifndef MUGY_DATA_DEV
#define MUGY_DATA_DEV

#include <cuComplex.h>  /* For complex data types. */
#include "mh_userFLAGS.h"
#include "mh_macros.h"

// Number of dimensions in the code.
#define nDim 3

#if USE_SINGLE_PRECISION > 0
typedef float real;
typedef cuComplex cufourier;
#else
typedef double real;
typedef cuDoubleComplex cufourier;
#endif

// Define our own int in case we wish to change to long.
typedef int mint;

// Moment indices.
#define denIdx 0
#define tempIdx 1

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Copy real-space moments between host and device.
//void memcpy_real_dev(real *dest, real *src, mint dofs, enum memcpy_dir_dev dir);
void memcpy_real_dev(real *dest, real *src, mint dofs, enum cudaMemcpyKind dir);

// Copy fourier-space moments between host and device.
//void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum memcpy_dir_dev dir);
void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum cudaMemcpyKind dir);

#endif
