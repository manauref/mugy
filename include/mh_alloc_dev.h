/* mugy: mh_alloc_dev.h
 *
 * Functions used to allocate arrays on device (GPU).
 *
 */
#pragma once

#include "mh_macros.h"

real* mugy_alloc_real_dev(mint numElements);
void* mugy_alloc_fourier_dev(mint numElements);
void* mugy_alloc_dev(mint numElements, size_t elemsz);
void mugy_free_dev(void *arr);
