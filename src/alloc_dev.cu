/* mugy: alloc_dev.cu
 *
 * Functions used to allocate arrays on device (GPU).
 *
 */

extern "C" {
#include "mh_alloc_dev.h"
#include "mh_fourier_dev.h"
#include "mh_utilities_dev.h"
}

real* mugy_alloc_real_dev(int numElements) {
  real *out_p;
  checkCudaErrors(cudaMalloc(&out_p, numElements*sizeof(real)));
  return out_p;
}
// MF 2023/03/29: Returning a void* because C doesn't know cuCumplex/mugy_cufourier_t the type.
void* mugy_alloc_fourier_dev(mint numElements) {
  mugy_cufourier_t *out_p;
  checkCudaErrors(cudaMalloc(&out_p, numElements*sizeof(mugy_cufourier_t)));
  return out_p;
}
void* mugy_alloc_dev(mint numElements, size_t elemsz) {
  void *out_p;
  checkCudaErrors(cudaMalloc(&out_p, numElements*elemsz));
  return out_p;
}

void mugy_free_dev(void *arr) {
  checkCudaErrors(cudaFree(arr));  // Free device memory.
}
