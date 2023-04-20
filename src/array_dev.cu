/* mugy: array_dev.cu
 *  
 * Device methods for mugy_arrays.
 *
 */

extern "C" {
#include "mh_array.h"
#include "mh_utilities.h"
#include "mh_array_dev.h"
#include "mh_data_dev.h"
}
#include "mh_fourier_dev.h"
#include "mh_utilities_dev.h"
#include <math.h>

// Starting linear index for each thread.
#define LINIDX0 (blockIdx.x*blockDim.x+threadIdx.x)

void mugy_array_zero_dev(struct mugy_array *arr) {
  // Set all elements in the array to zero.
  checkCudaErrors(cudaMemset(arr->dev, 0, arr->nelemsz));
}

__global__ void mugy_array_set_real_cu(real* out, mint nelem, real a, const real* inp) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = a*inp[linc];
}
__global__ void mugy_array_set_fourier_cu(mugy_cufourier_t* out, mint nelem, real a, const mugy_cufourier_t* inp) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCmul(mugy_make_cuComplex(a,0.), inp[linc]);
}
void mugy_array_scale_dev(struct mugy_array *arr, real fac) {
  // Scale mugy_array data by a constant.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(arr->nelem, nthreads);
  if (arr->type == MUGY_REAL)
    mugy_array_set_real_cu<<<nblocks, nthreads>>>((real *)arr->dev, arr->nelem, fac, (const real *)arr->dev);
  else if (arr->type == MUGY_FOURIER)
    mugy_array_set_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)arr->dev, arr->nelem, fac, (const mugy_cufourier_t *)arr->dev);
}

__global__ void mugy_array_isfinite_real_cu(bool *out, real* in, mint nelem) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x) {
    out[0] = isfinite(in[linc]);
    if (!out[0]) break;
  }
}
__global__ void mugy_array_isfinite_fourier_cu(bool *out, mugy_cufourier_t* in, mint nelem) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x) {
    out[0] = isfinite(mugy_cuCreal(in[linc])) && isfinite(mugy_cuCimag(in[linc]));;
    if (!out[0]) break;
  }
}
bool mugy_array_isfinite_dev(struct mugy_array *arr) {
  // Check that none of the elements are inf or NaN.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(arr->nelem, nthreads);
  if (arr->type == MUGY_REAL)
    mugy_array_isfinite_real_cu<<<nblocks, nthreads>>>((bool *)arr->bool_dev, (real *)arr->dev, arr->nelem);
  else if (arr->type == MUGY_FOURIER)
    mugy_array_isfinite_fourier_cu<<<nblocks, nthreads>>>((bool *)arr->bool_dev, (mugy_cufourier_t *)arr->dev, arr->nelem);

  mugy_memcpy_dev(arr->bool_ho, arr->bool_dev, sizeof(bool), MUGY_DEVICE2HOST);
  return ((bool*)arr->bool_ho)[0];
}
