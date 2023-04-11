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
#include "mh_array.h"
#include <stddef.h>

// Wrappers to basic functions that allocate plain arrays on host.
//   numElements: number of elements in the array.
mint* alloc_mintArray_ho(mint numElements);
char* alloc_charArray_ho(mint numElements);
real* alloc_realArray_ho(mint numElements);
void* alloc_fourierArray_ho(mint numElements);
void *mugy_alloc_ho(mint numElements, size_t elemsz);
void *mugy_alloc(mint numElements, size_t elemsz, enum resource_mem res);

void mugy_free_ho(void *arr);
void mugy_free(void *arr, enum resource_mem res);
#endif
