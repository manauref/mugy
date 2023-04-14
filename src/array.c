/* mugy: array.c
 *  
 * Methods to allocate, use and free mugy_arrays.
 *
 */

#include "mh_array.h"
#include "mh_array_dev.h"
#include "mh_alloc.h"
#include "mh_fourier_ho.h"
#include "mh_data.h"
#include <stdlib.h>  // for malloc.

// Allocate array on host, device, or both.
struct mugy_array *mugy_array_alloc(enum mugy_datatype type, mint numElements, enum resource_mem res) {

  struct mugy_array *arr = (struct mugy_array *) malloc(sizeof(struct mugy_array));
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

  return arr;
}

// Copy arrays betwen host and device, or within host or device.
void *mugy_array_copy(struct mugy_array *aout, struct mugy_array *ain, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir == host2device)
    return mugy_memcpy(aout->dev, ain->ho, ain->nelemsz, dir);
  else if (dir == device2host)
    return mugy_memcpy(aout->ho, ain->dev, ain->nelemsz, dir);
  else if (dir == device2device)
    return mugy_memcpy(aout->dev, ain->dev, ain->nelemsz, dir);
#endif

  return mugy_memcpy(aout->ho, ain->ho, ain->nelemsz, host2host);
}

// Scale array by a constant 'fac'.
void mugy_array_scale(struct mugy_array *arr, real fac, enum resource_comp res) {
#ifdef USE_GPU
  if (res == deviceComp)
    return mugy_array_scale_dev(arr, fac);
#endif

  real *ap = arr->ho;
  for (mint linIdx=0; linIdx<arr->nelem; linIdx++) {
    ap[0] *= fac;  ap++;
  }
}

// Free memory associated with arrays on host, device or both.
void mugy_array_free(struct mugy_array *arr, enum resource_mem res) {
  if ((res == hostMem) || (res == hostAndDeviceMem))
    mugy_free(arr->ho, hostMem);  // Free host memory.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    mugy_free(arr->dev, deviceMem);  // Free device memory.
}
