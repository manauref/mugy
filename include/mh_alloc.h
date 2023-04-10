/* mugy: mh_alloc.h
 *
 * Functions used to allocate arrays.
 *
 */

#ifndef MUGY_ALLOC
#define MUGY_ALLOC

#include "mh_data.h"
#include "mh_grid.h"
#include "mh_population.h"
#include "mh_utilities.h"

// Wrappers to basic functions that allocate plain arrays on host.
//   numElements: number of elements in the array.
mint* alloc_mintArray_ho(mint numElements);
char* alloc_charArray_ho(mint numElements);
real* alloc_realArray_ho(mint numElements);
fourier* alloc_fourierArray_ho(mint numElements);

// Functions that allocate real/Fourier arrays on host, device or both.
//   arr: pointer to array to be allocated (really a struct mugy_with ho/dev pointers).
//   numElements: number of elements in the array.
//   res: indicates which resource to use (host, device, both).
void alloc_realArray(struct mugy_realArray *arr, mint numElements, enum resource_mem res);
void alloc_fourierArray(struct mugy_fourierArray *arrk, mint numElements, enum resource_mem res);

/* Allocate real-space moment vectors, on host and/or device.
     mom: struct mugy_holding the vector of moments.
     grid: grid on which to allocate the vector of moments.
     pop: population struct mugy_containing the number of species and moments.
     res: resource on which to allocate (host, device or both). */
void alloc_realMoments(struct mugy_realArray *mom, const struct mugy_realGrid grid, const struct mugy_population pop, enum resource_mem res);

/* Allocate Fourier-space moment vectors, on host and/or device.
     momk: struct mugy_holding the vector of moments.
     grid: grid on which to allocate the vector of moments.
     pop: population struct mugy_containing the number of species and moments.
     res: resource on which to allocate (host, device or both). */
void alloc_fourierMoments(struct mugy_fourierArray *momk, const struct mugy_fourierGrid grid, const struct mugy_population pop, enum resource_mem res);

/* Functions to free memory associated used for moment vector. 
     arr: array to free.
     res: resource where memory needs to be freed (host, device or both). */
void free_realArray(struct mugy_realArray *arr, enum resource_mem res);
void free_fourierArray(struct mugy_fourierArray *arrk, enum resource_mem res);

#endif
