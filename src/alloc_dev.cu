/* mugy: alloc_dev.cu

   Functions used to allocate arrays on device (GPU).
*/

extern "C" {
#include "mh_utilities_dev.h"
#include "mh_data_dev.h"
#include "mh_alloc_dev.h"
}

real* alloc_realArray_dev(int numElements) {
  real *out_p;
  checkCudaErrors(cudaMalloc(&out_p, numElements*sizeof(real)));
  return out_p;
}
// MF 2023/03/29: Returning a void* because C doesn't know cuCumplex/cufourier the type.
void* alloc_fourierArray_dev(const mint numElements) {
  cufourier *out_p;
  checkCudaErrors(cudaMalloc(&out_p, numElements*sizeof(cufourier)));
  return out_p;
}

void free_realArray_dev(real *arr) {
  checkCudaErrors(cudaFree(arr));  // Free device memory.
}
void free_fourierArray_dev(void *arrk) {
  checkCudaErrors(cudaFree((cufourier*) arrk));  // Free device memory.
}