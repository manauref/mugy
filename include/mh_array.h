/* mugy: array.c
 *  
 * Methods to allocate, use and free mugy_arrays.
 *
 */
#ifndef MUGY_ARRAY
#define MUGY_ARRAY

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

// Functions that allocate real/Fourier arrays on host, device or both.
//   arr: pointer to array to be allocated (really a struct mugy_with ho/dev pointers).
//   numElements: number of elements in the array.
//   res: indicates which resource to use (host, device, both).
void mugy_array_alloc(struct mugy_array *arr, enum mugy_datatype type, mint numElements, enum resource_mem res);

// Function to free memory associated used for mugy_array. 
//   arr: array to free.
//   res: resource where memory needs to be freed (host, device or both).
void mugy_array_free(struct mugy_array *arr, enum resource_mem res);

// Copy real/fourier array between host and device.
void mugy_array_hodevXfer(struct mugy_array *arr, enum memcpy_dir_dev dir);

// Scale an array by a factor 'fac'.
void mugy_array_scale(struct mugy_array *arr, real fac, enum resource_comp res);

#endif
