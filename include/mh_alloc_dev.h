/* mugy: mh_alloc_dev.h

   Functions used to allocate arrays on device (GPU).
*/

#ifndef MUGY_ALLOC_DEV
#define MUGY_ALLOC_DEV

//#include "mh_alloc.h"
#include "mh_data_dev.h"

real* alloc_realArray_dev(mint numElements);
void* alloc_fourierArray_dev(mint numElements);

void free_realArray_dev(real *arr);
void free_fourierArray_dev(void *arrk);

#endif
