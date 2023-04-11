/* mugy: array_dev.cu
 *  
 * Device methods for mugy_arrays.
 *
 */

extern "C" {
#include "mh_array.h"
#include "mh_utilities.h"
#include "mh_array_dev.h"
}
#include "mh_fourier_dev.h"

// Starting linear index for each thread.
#define LINIDX0 (blockIdx.x*blockDim.x+threadIdx.x)

__global__ void
mugy_array_set_real_cu(real* out, mint nelem, real a, const real* inp)
{
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = a*inp[linc];
}

__global__ void
mugy_array_set_fourier_cu(cufourier* out, mint nelem, real a, const cufourier* inp)
{
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCmul(mugy_make_cuDoubleComplex(a,0.),inp[linc]);
}

// Scale mugy_array data by a constant.
void mugy_array_scale_dev(struct mugy_array *arr, real fac) {
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(arr->nelem, nthreads);
  if (arr->type == real_enum)
    mugy_array_set_real_cu<<<nblocks, nthreads>>>((real *)arr->dev, arr->nelem, fac, (const real *)arr->dev);
  else if (arr->type == fourier_enum)
    mugy_array_set_fourier_cu<<<nblocks, nthreads>>>((cufourier *)arr->dev, arr->nelem, fac, (const cufourier *)arr->dev);
}
