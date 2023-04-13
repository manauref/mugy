/* mugy: mh_alloc_dev.h
 *
 * Functions used to allocate arrays on device (GPU).
 *
 */
#pragma once

#include "mh_macros.h"

real* alloc_realArray_dev(mint numElements);
void* alloc_fourierArray_dev(mint numElements);
void *mugy_alloc_dev(mint numElements, size_t elemsz);

void free_realArray_dev(real *arr);
void free_fourierArray_dev(void *arrk);
void mugy_free_dev(void *arr);
