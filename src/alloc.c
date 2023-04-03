/* mugy: alloc
   
   Functions used to allocate arrays.
*/
#include "mh_data.h"
#include "mh_utilities.h"
#include <stdlib.h>  // e.g. for calloc.
#include "mh_alloc.h"
#include "mh_alloc_dev.h"

// Wrappers to basic functions that allocate memory.
mint* alloc_mintArray_ho(mint numElements) {
  mint *out_p;
  out_p = (mint *) calloc(numElements, sizeof(mint));
  if (out_p == NULL) abortSimulation(" alloc_mintArray: calloc failed! Terminating...\n");
  return out_p;
}
char* alloc_charArray_ho(mint numElements) {
  char *out_p;
  out_p = (char *) calloc(numElements, sizeof(char));
  if (out_p == NULL) abortSimulation(" alloc_charArray: calloc failed! Terminating...\n");
  return out_p;
}
real* alloc_realArray_ho(mint numElements) {
  real *out_p;
  out_p = (real *) calloc(numElements, sizeof(real));
  if (out_p == NULL)
    abortSimulation(" alloc_realArray: calloc failed! Terminating...\n");
  return out_p;
}
fourier *alloc_fourierArray_ho(mint numElements) {
  fourier *out_p;
  out_p = (fourier *) calloc(numElements, sizeof(fourier));
  if (out_p == NULL)
    abortSimulation(" alloc_fourierArray: calloc failed! Terminating...\n");
  return out_p;
}

// Functions that allocate arrays on host, device, or both.
void alloc_realArray(struct realArray *arr, mint numElements, enum resource_mem res) {
  arr->nelem = numElements;

  if ((res == hostMem) || (res == hostAndDeviceMem))
    arr->ho = alloc_realArray_ho(arr->nelem);  // Allocate on host.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    arr->dev = alloc_realArray_dev(arr->nelem);  // Allocate on device.
}
void alloc_fourierArray(struct fourierArray *arrk, mint numElements, enum resource_mem res) {
  arrk->nelem = numElements;

  if ((res == hostMem) || (res == hostAndDeviceMem))
    arrk->ho = alloc_fourierArray_ho(arrk->nelem);  // Allocate on host.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    arrk->dev = alloc_fourierArray_dev(arrk->nelem);  // Allocate on device.
}

// Functions that allocate moment vectors.
void alloc_realMoments(struct realArray *mom, const struct realGrid grid, const struct population pop, enum resource_mem res) {
  mint nelem = pop.numMomentsTot*prod_mint(grid.Nx,nDim);
  alloc_realArray(mom, nelem, res);
}
void alloc_fourierMoments(struct fourierArray *momk, const struct fourierGrid grid, const struct population pop, enum resource_mem res) {
  mint nelem = pop.numMomentsTot*prod_mint(grid.Nekx,nDim);
  alloc_fourierArray(momk, nelem, res);
}

// Functions to free memory associated with arrays on host, device or both.
void free_realArray(struct realArray *arr, enum resource_mem res) {
  if ((res == hostMem) || (res == hostAndDeviceMem))
    free(arr->ho);  // Free host memory.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    free_realArray_dev(arr->dev);  // Free device memory.
}
void free_fourierArray(struct fourierArray *arrk, enum resource_mem res) {
  if ((res == hostMem) || (res == hostAndDeviceMem))
    free(arrk->ho);  // Free host memory.

  if ((res == deviceMem) || (res == hostAndDeviceMem))
    free_fourierArray_dev(arrk->dev);  // Free device memory.
}
