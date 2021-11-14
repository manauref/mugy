/* mugy: alloc_mugy
   
   Functions used to allocate arrays.
*/
#include "alloc_mugy.h"

// Wrappers to basic functions that allocate memory.
int* alloc_intArray(const int numElements) {
  int *out_p;
  out_p = (int *) calloc(numElements, sizeof(int));
  if (out_p == NULL) abortSimulation(" alloc_intArray: calloc failed! Terminating...\n");
  return out_p;
}
char* alloc_charArray(const int numElements) {
  char *out_p;
  out_p = (char *) calloc(numElements, sizeof(char));
  if (out_p == NULL) abortSimulation(" alloc_charArray: calloc failed! Terminating...\n");
  return out_p;
}
real* alloc_realArray(const int numElements) {
  real *out_p;
  out_p = (real *) calloc(numElements, sizeof(real));
  if (out_p == NULL)
    abortSimulation(" alloc_realArray: calloc failed! Terminating...\n");
  return out_p;
}
fourier *alloc_fourierArray(const int numElements) {
  fourier *out_p;
  out_p = (fourier *) calloc(numElements, sizeof(fourier));
  if (out_p == NULL)
    abortSimulation(" alloc_fourierArray: calloc failed! Terminating...\n");
  return out_p;
}

// Functions that allocate moment vectors, on host and/or device.
void alloc_realMoments(const struct realGrid grid, const struct speciesParameters spec, const resource res, struct realMoments *mom) {
  if ((res == hostOnly) || (res == hostAndDevice))
    mom->ho = alloc_realArray(spec.numSpecies*spec.numMoments*prod_int(grid.Nx,nDim));  // Allocate on host.
}

void alloc_fourierMoments(const struct fourierGrid grid, const struct speciesParameters spec, const resource res, struct fourierMoments *momk) {
  if ((res == hostOnly) || (res == hostAndDevice))
    momk->ho = alloc_fourierArray(spec.numSpecies*spec.numMoments*prod_int(grid.Nekx,nDim));  // Allocate on host.
}

// Functions to free memory associated used for moment vector. 
void free_realMoments(struct realMoments *mom, const resource res) {
  if ((res == hostOnly) || (res == hostAndDevice))
    free(mom->ho);  // Free host memory.
}
void free_fourierMoments(struct fourierMoments *momk, const resource res) {
  if ((res == hostOnly) || (res == hostAndDevice))
    free(momk->ho);  // Free host memory.
}
