/* mugy: data_dev.cu

   Functions to manipulate data on the device. For example:
     - Copy between host and device.
*/

extern "C" {
#include "mh_utilities_dev.h"
#include "mh_data_dev.h"
}

// MF 2023/03/30: for some reason I can't use memcpy_dir_dev here.
//void memcpy_realMoments_dev(real *dest, real *src, mint dofs, enum memcpy_dir_dev dir) {
void memcpy_real_dev(real *dest, real *src, mint numElements, enum cudaMemcpyKind dir) {
  checkCudaErrors(cudaMemcpy(dest, src, numElements*sizeof(real), dir));
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/cufourier the type.
//void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum memcpy_dir_dev dir) {
void memcpy_fourier_dev(void *dest, void *src, mint numElements, enum cudaMemcpyKind dir) {
  checkCudaErrors(cudaMemcpy(dest, src, numElements*sizeof(cufourier), dir));
}
