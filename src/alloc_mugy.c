/* mugy
   
   Functions used to allocate arrays.
*/
#include "alloc_mugy.h"

// Wrappers to basic functions that allocate memory.
real* alloc_realArray(const int numElements) {
  real *pOut;
  pOut = (real *) calloc(numElements, sizeof(real));
  if (pOut == NULL) {
    abortSimulation(" alloc_realArray: calloc failed! Terminating...\n");
  }
  return pOut;
}
fourier *alloc_fourierArray(const int numElements) {
  fourier *pOut;
  pOut = (fourier *) calloc(numElements, sizeof(fourier));
  if (pOut == NULL) {
    abortSimulation(" alloc_fourierArray: calloc failed! Terminating...\n");
  }
  return pOut;
}

// Functions that allocate moment vectors, on host and/or device.
void alloc_realMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct realMoments *mom) {
  if ((res == hostOnly) || (res == hostAndDevice)) {
    mom->ho = alloc_realArray(ns*nmom*prod(grid.NxG));  // Allocate on host.
  }
}

void alloc_fourierMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct fourierMoments *momk) {
  if ((res == hostOnly) || (res == hostAndDevice)) {
    momk->ho = alloc_fourierArray(ns*nmom*prod(grid.NkxG));  // Allocate on host.
  }
}

void alloc_aliasedRealMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct realMoments *mom) {
  if ((res == hostOnly) || (res == hostAndDevice)) {
    mom->ho = alloc_realArray(ns*nmom*prod(grid.NxaG));  // Allocate on host.
  }
}

void alloc_aliasedFourierMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct fourierMoments *momk) {
  if ((res == hostOnly) || (res == hostAndDevice)) {
    momk->ho = alloc_fourierArray(ns*nmom*prod(grid.NkxaG));  // Allocate on host.
  }
}

// Functions to free memory associated used for moment vector. 
void free_realMoments(struct realMoments *mom, const resource res) {
  if ((res == hostOnly) || (res == hostAndDevice)) {
    free(mom->ho);  // Free host memory.
  }
}
void free_fourierMoments(struct fourierMoments *momk, const resource res) {
  if ((res == hostOnly) || (res == hostAndDevice)) {
    free(momk->ho);  // Free host memory.
  }
}
