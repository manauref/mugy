/* mugy: data_dev.cu

   Functions to manipulate data on the device. For example:
     - Copy between host and device.
*/

extern "C" {
#include "mh_utilities_dev.h"
#include "mh_data_dev.h"
}

//void memcpy_realMoments_dev(real *dest, real *src, mint dofs, enum memcpy_dir_dev dir) {
void memcpy_real_dev(real *dest, real *src, mint dofs, enum cudaMemcpyKind dir) {
  checkCudaErrors(cudaMemcpy(dest, src, dofs*sizeof(real), dir));
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/cufourier the type.
//void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum memcpy_dir_dev dir) {
void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum cudaMemcpyKind dir) {
  checkCudaErrors(cudaMemcpy(dest, src, dofs*sizeof(cufourier), dir));
}
