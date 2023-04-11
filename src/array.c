/* mugy: array.c
 *  
 * Methods to allocate, use and free mugy_arrays.
 *
 */

#include "mh_array.h"
#include "mh_alloc.h"
#include "mh_fourier_ho.h"
#include "mh_data.h"

// Allocate array on host, device, or both.
void mugy_array_alloc(struct mugy_array *arr, enum mugy_datatype type, mint numElements, enum resource_mem res) {
  arr->type  = type;
  arr->nelem = numElements;

  if (arr->type == mint_enum)         arr->elemsz = sizeof(mint);
  else if (arr->type == real_enum)    arr->elemsz = sizeof(real);
  else if (arr->type == fourier_enum) arr->elemsz = sizeof(fourier);
  arr->nelemsz = arr->nelem*arr->elemsz;

  if ((res == hostMem) || (res == hostAndDeviceMem))
    arr->ho = mugy_alloc(arr->nelem, arr->elemsz, hostMem);  // Allocate on host.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    arr->dev = mugy_alloc(arr->nelem, arr->elemsz, deviceMem);  // Allocate on device.
}

// Free memory associated with arrays on host, device or both.
void mugy_array_free(struct mugy_array *arr, enum resource_mem res) {
  if ((res == hostMem) || (res == hostAndDeviceMem))
    mugy_free(arr->ho, hostMem);  // Free host memory.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    mugy_free(arr->dev, deviceMem);  // Free device memory.
}

// Copy real(Fourier)Arrays betwen host and device.
void mugy_array_hodevXfer(struct mugy_array *arr, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir == host2device)
    mugy_memcpy(arr->dev, arr->ho, arr->nelemsz, dir);
  else if (dir == device2host)
    mugy_memcpy(arr->ho, arr->dev, arr->nelemsz, dir);
#endif
}

// Functions that scale real(Fourier)Arrays by a constant.
void mugy_array_scale(struct mugy_array *arr, real fac, enum resource_comp res) {
//#ifdef USE_GPU
//  if (res == deviceComp)
//    return mugy_array_scale_dev(arr, fac, res);
//#endif

  real *fk = arr->ho;
  for (mint linIdx=0; linIdx<arr->nelem; linIdx++) {
    fk[0] *= fac;  fk++;
  }
}
