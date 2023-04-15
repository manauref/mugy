/* mugy: mh_array.h
 *  
 * Methods to allocate, use and free mugy_arrays.
 *
 */
#pragma once

#include "mh_macros.h"
#include <stddef.h>

// Structure storing an array on host, device, or both.
struct mugy_array {
  enum mugy_datatype type;  // Type of the data (e.g. real, fourier).
  mint nelem;  // Number of elements allocated.
  size_t elemsz;  // Size of each element (in bytes, as returned by sizeof).
  size_t nelemsz; // nelem*elemsz;
  void *ho;    // Pointer to host memory.
  void *dev;   // Pointer to device memory.
};

// Get pointer to the element in the array with linear index 'linIdx'.
static inline void* mugy_array_get(struct mugy_array* arr, mint linIdx) {
  return ((char *)arr->ho) + linIdx*arr->elemsz;
};

// Functions that allocate real/Fourier arrays on host, device or both.
//   arr: pointer to array to be allocated (really a struct mugy_with ho/dev pointers).
//   numElements: number of elements in the array.
//   res: indicates which resource to use (host, device, both).
struct mugy_array *mugy_array_alloc(enum mugy_datatype type, mint numElements, enum resource_mem res);

// Function to free memory associated used for mugy_array. 
//   arr: array to free.
//   res: resource where memory needs to be freed (host, device or both).
void mugy_array_free(struct mugy_array *arr, enum resource_mem res);

// Copy arrays betwen host and device, or within host or device.
void *mugy_array_copy(struct mugy_array *aout, struct mugy_array *ain, enum memcpy_dir_dev dir);

// Scale an array by a factor 'fac'.
void mugy_array_scale(struct mugy_array *arr, real fac, enum resource_comp res);
