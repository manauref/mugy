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
#include <string.h>  // for memset.
#include <math.h>  // for isfinite.
#include <assert.h> // for assert.

struct mugy_array *mugy_array_alloc(enum mugy_data_types type, mint numElements, enum mugy_resource_mem res) {
  // Allocate array on host, device, or both.

  struct mugy_array *arr = (struct mugy_array *) malloc(sizeof(struct mugy_array));
  arr->type  = type;
  arr->nelem = numElements;

  if (arr->type == MUGY_MINT)         arr->elemsz = sizeof(mint);
  else if (arr->type == MUGY_REAL)    arr->elemsz = sizeof(real);
  else if (arr->type == MUGY_FOURIER) arr->elemsz = sizeof(fourier);
  arr->nelemsz = arr->nelem*arr->elemsz;

  if ((res == MUGY_HOST_MEM) || (res == MUGY_HOSTDEVICE_MEM)) {
    // Allocate on host.
    arr->ho = mugy_alloc(arr->nelem, arr->elemsz, MUGY_HOST_MEM);
    arr->binflag_ho = mugy_alloc(1, sizeof(mint), MUGY_HOST_MEM);
  }

  if ((res == MUGY_DEVICE_MEM) || (res == MUGY_HOSTDEVICE_MEM)) {
    // Allocate on device.
    arr->dev = mugy_alloc(arr->nelem, arr->elemsz, MUGY_DEVICE_MEM);
    arr->binflag_dev = mugy_alloc(1, sizeof(mint), MUGY_DEVICE_MEM);
  }

  return arr;
}

void mugy_array_zero(struct mugy_array *arr, enum mugy_resource_mem res) {
  // Set all elements in the array to zero.
#ifdef USE_GPU
  if (res == MUGY_DEVICE_MEM) {
    mugy_array_zero_dev(arr);
    return;
  } else if (res == MUGY_HOSTDEVICE_MEM) {
    mugy_array_zero_dev(arr);
  }
#endif
  memset(arr->ho, 0, arr->nelemsz);
}

void *mugy_array_copy(struct mugy_array *aout, struct mugy_array *ain, enum mugy_memcpy_dir dir) {
  // Copy arrays betwen host and device, or within host or device.
#ifdef USE_GPU
  if (dir == MUGY_HOST2DEVICE)
    return mugy_memcpy(aout->dev, ain->ho, ain->nelemsz, dir);
  else if (dir == MUGY_DEVICE2HOST)
    return mugy_memcpy(aout->ho, ain->dev, ain->nelemsz, dir);
  else if (dir == MUGY_DEVICE2DEVICE)
    return mugy_memcpy(aout->dev, ain->dev, ain->nelemsz, dir);
#endif

  return mugy_memcpy(aout->ho, ain->ho, ain->nelemsz, MUGY_HOST2HOST);
}

void mugy_array_increment(struct mugy_array *out, struct mugy_array *x, enum mugy_resource_calc res) {
  // Increment array out by array x.
  assert(out->nelem == x->nelem);  assert(out->type == x->type);

#ifdef USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_array_increment_dev(out, x);
#endif

  real *outp = out->ho, *xp = x->ho;
  for (mint linIdx=0; linIdx<out->nelem; linIdx++) {
    outp[0] += xp[0];
    outp++;  xp++;
  }
}

void mugy_array_ax_assign(struct mugy_array *out, real a, struct mugy_array *x, enum mugy_resource_calc res) {
  // Assign array out with a*x.
  assert(out->nelem == x->nelem);  assert(out->type == x->type);

#ifdef USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_array_ax_assign_dev(out, a, x);
#endif

  real *outp = out->ho, *xp = x->ho;
  for (mint linIdx=0; linIdx<out->nelem; linIdx++) {
    outp[0] = a*xp[0];
    outp++;  xp++;
  }
}

void mugy_array_axpy_assign(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y, enum mugy_resource_calc res) {
  // Assign array out with a*x+y.
  assert((out->nelem == x->nelem) && (out->nelem == y->nelem));
  assert((out->type == x->type) && (out->type == y->type));

#ifdef USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_array_axpy_assign_dev(out, a, x, y);
#endif

  real *outp = out->ho, *xp = x->ho, *yp = y->ho;
  for (mint linIdx=0; linIdx<out->nelem; linIdx++) {
    outp[0] = a*xp[0]+yp[0];
    outp++;  xp++;  yp++;
  }
}

void mugy_array_ax_increment(struct mugy_array *out, real a, struct mugy_array *x, enum mugy_resource_calc res) {
  // Increment array out by a*x.
  mugy_array_axpy_assign(out, a, x, out, res);
}

void mugy_array_axpy_increment(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y, enum mugy_resource_calc res) {
  // Increment array out by a*x+y.
  assert((out->nelem == x->nelem) && (out->nelem == y->nelem));
  assert((out->type == x->type) && (out->type == y->type));

#ifdef USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_array_axpy_increment_dev(out, a, x, y);
#endif

  real *outp = out->ho, *xp = x->ho, *yp = y->ho;
  for (mint linIdx=0; linIdx<out->nelem; linIdx++) {
    outp[0] += a*xp[0]+yp[0];
    outp++;  xp++;  yp++;
  }
}

void mugy_array_scale(struct mugy_array *arr, real fac, enum mugy_resource_calc res) {
  // Scale array by a constant 'fac'.
#ifdef USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_array_scale_dev(arr, fac);
#endif

  real *ap = arr->ho;
  for (mint linIdx=0; linIdx<arr->nelem; linIdx++) {
    ap[0] *= fac;  ap++;
  }
}

bool mugy_array_isfinite(struct mugy_array *arr, enum mugy_resource_calc res) {
  // Check that none of the elements are inf or NaN.
#ifdef USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_array_isfinite_dev(arr);
#endif

  real *ap = arr->ho;
  bool isvalid = true;
  for (mint linIdx=0; linIdx<arr->nelem; linIdx++) {
    isvalid = isfinite(ap[0]);
    if (!isvalid) break;
    ap++;
  }
  return isvalid;
}

// Free memory associated with arrays on host, device or both.
void mugy_array_free(struct mugy_array *arr, enum mugy_resource_mem res) {
  if ((res == MUGY_HOST_MEM) || (res == MUGY_HOSTDEVICE_MEM)) {
    // Free host memory.
    mugy_free(arr->ho, MUGY_HOST_MEM);
    mugy_free(arr->binflag_ho, MUGY_HOST_MEM);
  }

  if ((res == MUGY_DEVICE_MEM) || (res == MUGY_HOSTDEVICE_MEM)) {
    // Free device memory.
    mugy_free(arr->dev, MUGY_DEVICE_MEM);
    mugy_free(arr->binflag_dev, MUGY_DEVICE_MEM);
  }
}
