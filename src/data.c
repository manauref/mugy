/* mugy: data.c
   
   Parameters and fields used throughout mugy.
*/

#include "mh_data.h"
#include "mh_data_dev.h"
#include <string.h>

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

// Functions that copy memory between host and device.
void memcpy_mint(mint *dest, mint *src, mint numElements, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir != host2host)
    return memcpy_mint_dev(dest, src, numElements, dir);
#endif
  memcpy_mint_ho(dest, src, numElements);
}
void memcpy_real(real *dest, real *src, mint numElements, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir != host2host)
    return memcpy_real_dev(dest, src, numElements, dir);
#endif
  memcpy_real_ho(dest, src, numElements);
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/cufourier the type.
void memcpy_fourier(void *dest, void *src, mint numElements, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir != host2host)
    return memcpy_fourier_dev(dest, src, numElements, dir);
#endif
  memcpy_fourier_ho(dest, src, numElements);
}

// Functions that copy real(Fourier)Arrays betwen host and device.
void hodevXfer_realArray(struct mugy_realArray *arr, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir == host2device)
    memcpy_real_dev(arr->dev, arr->ho, arr->nelem, dir);
  else if (dir == device2host)
    memcpy_real_dev(arr->ho, arr->dev, arr->nelem, dir);
#endif
}
void hodevXfer_fourierArray(struct mugy_fourierArray *arr, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir == host2device)
    memcpy_fourier_dev(arr->dev, arr->ho, arr->nelem, dir);
  else if (dir == device2host)
    memcpy_fourier_dev(arr->ho, arr->dev, arr->nelem, dir);
#endif
}

// Functions that scale real(Fourier)Arrays by a constant.
void scale_realArray(struct mugy_realArray *arr, real fac, enum resource_comp res) {
#ifdef USE_GPU
//  if (res == deviceComp)
//    return scale_realArray_dev(arr, fac, res);
#endif

  real *fk = arr->ho;
  for (mint linIdx=0; linIdx<arr->nelem; linIdx++) {
    fk[0] *= fac;  fk++;
  }
}
void scale_fourierArray(struct mugy_fourierArray *arrk, real fac, enum resource_comp res) {
#ifdef USE_GPU
//  if (res == deviceComp)
//    return scale_fourierArray_dev(arr, fac, res);
#endif

  fourier *fk = arrk->ho;
  for (mint linIdx=0; linIdx<arrk->nelem; linIdx++) {
    fk[0] *= fac;  fk++;
  }
}
