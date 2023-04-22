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

__global__ void mugy_array_increment_real_cu(real *out, mint nelem, const real *inp) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] += inp[linc];
}
__global__ void mugy_array_increment_fourier_cu(mugy_cufourier_t *out, mint nelem, const mugy_cufourier_t *inp) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCadd(out[linc], inp[linc]);
}
void mugy_array_increment_dev(struct mugy_array *out, struct mugy_array *x) {
  // Increment array out by x.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(out->nelem, nthreads);
  if (out->type == MUGY_REAL)
    mugy_array_increment_real_cu<<<nblocks, nthreads>>>((real *)out->dev, out->nelem,
      (const real *)x->dev);
  else if (out->type == MUGY_FOURIER)
    mugy_array_increment_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)out->dev, out->nelem,
      (const mugy_cufourier_t *)x->dev);
}

__global__ void mugy_array_ax_assign_real_cu(real *out, mint nelem, real a, const real *inp) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = a*inp[linc];
}
__global__ void mugy_array_ax_assign_fourier_cu(mugy_cufourier_t *out, mint nelem, real a, const mugy_cufourier_t *inp) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCmul(mugy_make_cuComplex(a,0.), inp[linc]);
}
void mugy_array_ax_assign_dev(struct mugy_array *out, real a, struct mugy_array *x) {
  // Assign array out with a*x.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(out->nelem, nthreads);
  if (out->type == MUGY_REAL)
    mugy_array_ax_assign_real_cu<<<nblocks, nthreads>>>((real *)out->dev, out->nelem,
      a, (const real *)x->dev);
  else if (out->type == MUGY_FOURIER)
    mugy_array_ax_assign_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)out->dev, out->nelem,
      a, (const mugy_cufourier_t *)x->dev);
}

__global__ void mugy_array_axpy_assign_real_cu(real *out, mint nelem, real a, const real *x, const real *y) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = a*x[linc]+y[linc];
}
__global__ void mugy_array_axpy_assign_fourier_cu(mugy_cufourier_t *out, mint nelem, real a, const mugy_cufourier_t *x, const mugy_cufourier_t *y) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCadd( mugy_cuCmul(mugy_make_cuComplex(a,0.), x[linc]), y[linc]);
}
void mugy_array_axpy_assign_dev(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y) {
  // Assign array out with a*x+y.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(out->nelem, nthreads);
  if (out->type == MUGY_REAL)
    mugy_array_axpy_assign_real_cu<<<nblocks, nthreads>>>((real *)out->dev, out->nelem,
      a, (const real *)x->dev, (const real *)y->dev);
  else if (out->type == MUGY_FOURIER)
    mugy_array_axpy_assign_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)out->dev, out->nelem,
      a, (const mugy_cufourier_t *)x->dev, (const mugy_cufourier_t *)y->dev);
}

__global__ void mugy_array_axpy_increment_real_cu(real *out, mint nelem, real a, const real *x, const real *y) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] += a*x[linc]+y[linc];
}
__global__ void mugy_array_axpy_increment_fourier_cu(mugy_cufourier_t *out, mint nelem, real a, const mugy_cufourier_t *x, const mugy_cufourier_t *y) {
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x)
    out[linc] = mugy_cuCadd( out[linc], mugy_cuCadd( mugy_cuCmul(mugy_make_cuComplex(a,0.), x[linc]), y[linc]) );
}
void mugy_array_axpy_increment_dev(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y) {
  // Increment array out with a*x+y.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(out->nelem, nthreads);
  if (out->type == MUGY_REAL)
    mugy_array_axpy_increment_real_cu<<<nblocks, nthreads>>>((real *)out->dev, out->nelem,
      a, (const real *)x->dev, (const real *)y->dev);
  else if (out->type == MUGY_FOURIER)
    mugy_array_axpy_increment_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)out->dev, out->nelem,
      a, (const mugy_cufourier_t *)x->dev, (const mugy_cufourier_t *)y->dev);
}

void mugy_array_scale_dev(struct mugy_array *arr, real fac) {
  // Scale mugy_array data by a constant.
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(arr->nelem, nthreads);
  if (arr->type == MUGY_REAL)
    mugy_array_ax_assign_real_cu<<<nblocks, nthreads>>>((real *)arr->dev, arr->nelem, fac, (const real *)arr->dev);
  else if (arr->type == MUGY_FOURIER)
    mugy_array_ax_assign_fourier_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)arr->dev, arr->nelem, fac, (const mugy_cufourier_t *)arr->dev);
}

__global__ void mugy_array_isnotfinite_real_cu(mint *out, real* in, mint nelem) {
  bool validval = true;
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x) {
    validval = isfinite(in[linc]);
    if (!validval) break;
  }
  if (!validval) atomicExch(out, 1);
}
__global__ void mugy_array_isnotfinite_fourier_cu(mint *out, mugy_cufourier_t* in, mint nelem) {
  bool validval = true;
  for (unsigned long linc = LINIDX0; linc < nelem; linc += blockDim.x*gridDim.x) {
    validval = isfinite(mugy_cuCreal(in[linc])) && isfinite(mugy_cuCimag(in[linc]));
    if (!validval) break;
  }
  if (!validval) atomicExch(out, 1);
}
bool mugy_array_isfinite_dev(struct mugy_array *arr) {
  // Check that none of the elements are inf or NaN.
  checkCudaErrors(cudaMemset(arr->binflag_dev, 0, sizeof(mint)));

  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(arr->nelem, nthreads);
  if (arr->type == MUGY_REAL)
    mugy_array_isnotfinite_real_cu<<<nblocks, nthreads>>>((mint *)arr->binflag_dev, (real *)arr->dev, arr->nelem);
  else if (arr->type == MUGY_FOURIER)
    mugy_array_isnotfinite_fourier_cu<<<nblocks, nthreads>>>((mint *)arr->binflag_dev, (mugy_cufourier_t *)arr->dev, arr->nelem);

  mugy_memcpy_dev(arr->binflag_ho, arr->binflag_dev, sizeof(mint), MUGY_DEVICE2HOST);
  return arr->binflag_ho[0]==0;
}
