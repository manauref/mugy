/* mugy: mh_array.h
 *  
 * Methods to allocate, use and free mugy_arrays.
 *
 */
#pragma once

#include "mh_macros.h"
#include <stddef.h>
#include <stdbool.h>

// Structure storing an array on host, device, or both.
struct mugy_array {
  enum mugy_data_types type;  // Type of the data (e.g. real, fourier).
  mint nelem;  // Number of elements allocated.
  size_t elemsz;  // Size of each element (in bytes, as returned by sizeof).
  size_t nelemsz; // nelem*elemsz;
  void *ho;    // Pointer to host memory.
  void *dev;   // Pointer to device memory.
  // Device memory to store the result of some reductions.
  mint *binflag_ho, *binflag_dev;
};

// Functions that allocate real/Fourier arrays on host, device or both.
//   arr: pointer to array to be allocated (really a struct mugy_with ho/dev pointers).
//   numElements: number of elements in the array.
//   res: indicates which resource to use (host, device, both).
struct mugy_array *mugy_array_alloc(enum mugy_data_types type, mint numElements, enum mugy_resource_mem res);

// Get pointer to the element in the array with linear index 'linIdx'.
static inline void* mugy_array_get(struct mugy_array* arr, mint linIdx) {
  return ((char *)arr->ho) + linIdx*arr->elemsz;
};

// Set all elements in the array to zero.
void mugy_array_zero(struct mugy_array *arr, enum mugy_resource_mem res);

// Function to free memory associated used for mugy_array. 
//   arr: array to free.
//   res: resource where memory needs to be freed (host, device or both).
void mugy_array_free(struct mugy_array *arr, enum mugy_resource_mem res);

// Copy arrays betwen host and device, or within host or device.
void *mugy_array_copy(struct mugy_array *aout, struct mugy_array *ain, enum mugy_memcpy_dir dir);

// Increment array out by array x.
void mugy_array_increment(struct mugy_array *out, struct mugy_array *x, enum mugy_resource_calc res);

// Assign array out with a*x.
void mugy_array_ax_assign(struct mugy_array *out, real a, struct mugy_array *x, enum mugy_resource_calc res);

// Assign array out with a*x+y.
void mugy_array_axpy_assign(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y, enum mugy_resource_calc res);

// Increment array out with a*x.
void mugy_array_ax_increment(struct mugy_array *out, real a, struct mugy_array *x, enum mugy_resource_calc res);

// Increment array out with a*x+y.
void mugy_array_axpy_increment(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y, enum mugy_resource_calc res);

// Scale an array by a factor 'fac'.
void mugy_array_scale(struct mugy_array *arr, real fac, enum mugy_resource_calc res);

// Check that none of the elements are inf or NaN.
bool mugy_array_isfinite(struct mugy_array *arr, enum mugy_resource_calc res);
