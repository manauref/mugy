/* mugy: data.c
 *  
 * Parameters and fields used throughout mugy.
 *
 */

#include "mh_data.h"
#include "mh_data_dev.h"
#include <string.h>  // For memcpy.
#include "mh_fourier_ho.h"

// Functions that copy memory on the host.
void memcpy_mint_ho(mint *dest, mint *src, mint numElements) {
  memcpy(dest, src, numElements*sizeof(mint));
}
void memcpy_real_ho(real *dest, real *src, mint numElements) {
  memcpy(dest, src, numElements*sizeof(real));
}
void memcpy_fourier_ho(fourier *dest, fourier *src, mint numElements) {
  memcpy(dest, src, numElements*sizeof(fourier));
}
void *mugy_memcpy_ho(void *dest, void *src, size_t sz) {
  return memcpy(dest, src, sz);
}

// Functions that copy memory between host and device.
void memcpy_mint(mint *dest, mint *src, mint numElements, enum mugy_memcpy_dir dir) {
#ifdef USE_GPU
  if (dir != MUGY_HOST2HOST)
    return memcpy_mint_dev(dest, src, numElements, dir);
#endif
  memcpy_mint_ho(dest, src, numElements);
}
void memcpy_real(real *dest, real *src, mint numElements, enum mugy_memcpy_dir dir) {
#ifdef USE_GPU
  if (dir != MUGY_HOST2HOST)
    return memcpy_real_dev(dest, src, numElements, dir);
#endif
  memcpy_real_ho(dest, src, numElements);
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/mugy_cufourier_t the type.
void memcpy_fourier(void *dest, void *src, mint numElements, enum mugy_memcpy_dir dir) {
#ifdef USE_GPU
  if (dir != MUGY_HOST2HOST)
    return memcpy_fourier_dev(dest, src, numElements, dir);
#endif
  memcpy_fourier_ho(dest, src, numElements);
}

void *mugy_memcpy(void *dest, void *src, size_t sz, enum mugy_memcpy_dir dir) {
#ifdef USE_GPU
  if (dir != MUGY_HOST2HOST)
    return mugy_memcpy_dev(dest, src, sz, dir);
#endif
  return mugy_memcpy_ho(dest, src, sz);
}
