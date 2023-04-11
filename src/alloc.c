/* mugy: alloc.c
 *
 * Functions used to allocate arrays.
 *
 */
#include "mh_alloc.h"
#include "mh_alloc_dev.h"
#include "mh_fourier_ho.h"
#include "mh_utilities.h"
#include <stdlib.h>  // e.g. for calloc.

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
  if (out_p == NULL) abortSimulation(" alloc_realArray: calloc failed! Terminating...\n");
  return out_p;
}
void *alloc_fourierArray_ho(mint numElements) {
  fourier *out_p;
  out_p = (fourier *) calloc(numElements, sizeof(fourier));
  if (out_p == NULL) abortSimulation(" alloc_fourierArray: calloc failed! Terminating...\n");
  return out_p;
}
void *mugy_alloc_ho(mint numElements, size_t elemsz) {
  void *out_p = calloc(numElements, elemsz);
  if (out_p == NULL) abortSimulation(" alloc_fourierArray: calloc failed! Terminating...\n");
  return out_p;
}
void *mugy_alloc(mint numElements, size_t elemsz, enum resource_mem res) {
  void *out_p;
  if (res == hostMem)
    out_p = mugy_alloc_ho(numElements, elemsz);  // Allocate on host.
  else if (res == deviceMem)
    out_p = mugy_alloc_dev(numElements, elemsz);  // Allocate on device.
  else
    abortSimulation(" mugy_alloc: invalid resource! Terminating...\n");
  return out_p;
}

void mugy_free_ho(void *arr) {
  free(arr);  // Deallocate on host.
}

void mugy_free(void *arr, enum resource_mem res) {
  if (res == hostMem)
    mugy_free_ho(arr);  // Deallocate on host.
  else if (res == deviceMem)
    mugy_free_dev(arr);  // Deallocate on device.
  else
    abortSimulation(" mugy_free: invalid resource! Terminating...\n");
}	
