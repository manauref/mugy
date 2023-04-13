/* mugy: mh_cuComplex_wrap.h
 *
 * Wrap cuComplex header so we can select the parts of it
 * we needed depending on precision.
 *
 */
#pragma once

#if (__NVCC__ && USE_GPU)

#include "mh_userFLAGS.h"

#include <cuComplex.h>  /* For complex data types. */

#ifdef USE_SINGLE_PRECISION


__host__ __device__ static __inline__ double mugy_cuCreal (cuFloatComplex x)
{
  return cuCrealf (x);
}

__host__ __device__ static __inline__ double mugy_cuCimag (cuFloatComplex x)
{
  return cuCimagf (x);
}

__host__ __device__ static __inline__ cuFloatComplex mugy_make_cuComplex
                                                           (double r, double i)
{
  return make_cuFloatComplex (r, i);
}

__host__ __device__ static __inline__ cuFloatComplex mugy_cuConj(cuFloatComplex x)
{
  return cuConjf(x);
}

__host__ __device__ static __inline__ cuFloatComplex mugy_cuCadd(cuFloatComplex x,
                                                             cuFloatComplex y)
{
  return cuCaddf(x, y);
}

__host__ __device__ static __inline__ cuFloatComplex mugy_cuCsub(cuFloatComplex x,
                                                             cuFloatComplex y)
{
  return cuCsubf(x, y);
}

__host__ __device__ static __inline__ cuFloatComplex mugy_cuCmul(cuFloatComplex x,
                                                             cuFloatComplex y)
{
  return cuCmulf(x, y);
}

__host__ __device__ static __inline__ cuFloatComplex mugy_cuCdiv(cuFloatComplex x,
                                                             cuFloatComplex y)
{
  return cuCdivf(x, y);
}

__host__ __device__ static __inline__ double mugy_cuCabs (cuFloatComplex x)
{
  return cuCabsf (x);
}

__host__ __device__ static __inline__  cuFloatComplex mugy_cuCfma( cuFloatComplex x, cuFloatComplex y, cuFloatComplex d)
{
  // Computes d+x*y.
  return cuCfmaf(x, y, d);
}


#else
// Double precision


__host__ __device__ static __inline__ double mugy_cuCreal (cuDoubleComplex x)
{
  return cuCreal (x);
}

__host__ __device__ static __inline__ double mugy_cuCimag (cuDoubleComplex x)
{
  return cuCimag (x);
}

__host__ __device__ static __inline__ cuDoubleComplex mugy_make_cuComplex
                                                           (double r, double i)
{
  return make_cuDoubleComplex (r, i);
}

__host__ __device__ static __inline__ cuDoubleComplex mugy_cuConj(cuDoubleComplex x)
{
  return cuConj(x);
}

__host__ __device__ static __inline__ cuDoubleComplex mugy_cuCadd(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
  return cuCadd(x, y);
}

__host__ __device__ static __inline__ cuDoubleComplex mugy_cuCsub(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
  return cuCsub(x, y);
}

__host__ __device__ static __inline__ cuDoubleComplex mugy_cuCmul(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
  return cuCmul(x, y);
}

__host__ __device__ static __inline__ cuDoubleComplex mugy_cuCdiv(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
  return cuCdiv(x, y);
}

__host__ __device__ static __inline__ double mugy_cuCabs (cuDoubleComplex x)
{
  return cuCabs (x);
}

__host__ __device__ static __inline__  cuDoubleComplex mugy_cuCfma( cuDoubleComplex x, cuDoubleComplex y, cuDoubleComplex d)
{
  // Computes d+x*y.
  return cuCfma(x, y, d);
}


#endif

#endif
