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
#include "mh_utilities_dev.h"

// Starting linear index for each thread.
#define LINIDX0 (blockIdx.x*blockDim.x+threadIdx.x)

void mugy_array_zero_dev(struct mugy_array *arr) {
  // Set all elements in the array to zero.
  checkCudaErrors(cudaMemset(arr->dev, 0, arr->nelemsz));
}

__global__ void
mugy_array_set_real_cu(real* out, mint nelem, real a, const real* inp)
{
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = a*inp[linc];
}

__global__ void
mugy_array_set_fourier_cu(mugy_cufourier_t* out, mint nelem, real a, const mugy_cufourier_t* inp)
{
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCmul(mugy_make_cuComplex(a,0.), inp[linc]);
}

// Scale mugy_array data by a constant.
void mugy_array_scale_dev(struct mugy_array *arr, real fac) {
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(arr->nelem, nthreads);
  if (arr->type == MUGY_REAL)
    mugy_array_set_real_cu<<<nblocks, nthreads>>>((real *)arr->dev, arr->nelem, fac, (const real *)arr->dev);
  else if (arr->type == MUGY_FOURIER)
    mugy_array_set_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)arr->dev, arr->nelem, fac, (const mugy_cufourier_t *)arr->dev);
}
