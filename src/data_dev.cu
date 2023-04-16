/* mugy: data_dev.cu
 *
 * Functions to manipulate data on the device. For example:
 *   - Copy between host and device.
 *
 */

extern "C" {
#include "mh_utilities_dev.h"
#include "mh_data_dev.h"
#include "mh_fourier_dev.h"
}

// MF 2023/03/30: for some reason we can't use memcpy_dir_dev directly cudaMemcpy.
// Maybe because of how nvcc handles an enum in C/C++/CUDA. Cast instead.
void memcpy_mint_dev(mint *dest, mint *src, mint numElements, enum memcpy_dir_dev dir) {
  checkCudaErrors(cudaMemcpy(dest, src, numElements*sizeof(mint), (cudaMemcpyKind) dir));
}
void memcpy_real_dev(real *dest, real *src, mint numElements, enum memcpy_dir_dev dir) {
  checkCudaErrors(cudaMemcpy(dest, src, numElements*sizeof(real), (cudaMemcpyKind) dir));
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/cufourier the type.
void memcpy_fourier_dev(void *dest, void *src, mint numElements, enum memcpy_dir_dev dir) {
  checkCudaErrors(cudaMemcpy(dest, src, numElements*sizeof(cufourier), (cudaMemcpyKind) dir));
}
void *mugy_memcpy_dev(void *dest, void *src, size_t sz, enum memcpy_dir_dev dir) {
  checkCudaErrors(cudaMemcpy(dest, src, sz, (cudaMemcpyKind) dir));
  return dest;
}
