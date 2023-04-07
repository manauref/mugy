/* mugy: mh_data_dev.h

   Data types and some macros used in mugy device (GPU) operations.
*/
#ifndef MUGY_DATA_DEV
#define MUGY_DATA_DEV

#include "mh_userFLAGS.h"
#include "mh_macros.h"

#if USE_GPU

#include <cuComplex.h>  /* For complex data types. */

#if USE_SINGLE_PRECISION
//typedef float real;
typedef cuComplex cufourier;
#else
//typedef double real;
typedef cuDoubleComplex cufourier;
#endif

#endif

// Copy mint/real/fourier data moments between host and device.
void memcpy_mint_dev(mint *dest, mint *src, mint dofs, enum memcpy_dir_dev dir);
void memcpy_real_dev(real *dest, real *src, mint dofs, enum memcpy_dir_dev dir);
void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum memcpy_dir_dev dir);

#endif
