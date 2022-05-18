/* mugy: alloc_mugy
   
   Functions used to allocate arrays.
*/
#include "mh_alloc.h"

// Wrappers to basic functions that allocate memory.
mint* alloc_mintArray(const mint numElements) {
  mint *out_p;
  out_p = (mint *) calloc(numElements, sizeof(mint));
  if (out_p == NULL) abortSimulation(" alloc_mintArray: calloc failed! Terminating...\n");
  return out_p;
}
char* alloc_charArray(const mint numElements) {
  char *out_p;
  out_p = (char *) calloc(numElements, sizeof(char));
  if (out_p == NULL) abortSimulation(" alloc_charArray: calloc failed! Terminating...\n");
  return out_p;
}
real* alloc_realArray(const mint numElements) {
  real *out_p;
  out_p = (real *) calloc(numElements, sizeof(real));
  if (out_p == NULL)
    abortSimulation(" alloc_realArray: calloc failed! Terminating...\n");
  return out_p;
}
fourier *alloc_fourierArray(const mint numElements) {
  fourier *out_p;
  out_p = (fourier *) calloc(numElements, sizeof(fourier));
  if (out_p == NULL)
    abortSimulation(" alloc_fourierArray: calloc failed! Terminating...\n");
  return out_p;
}

// Functions that allocate moment vectors, on host and/or device.
void alloc_realMoments(const struct realGrid grid, const struct population pop, const resource res, struct realMoments *mom) {
  mint numMomentsTot = 1;
  for (mint s=0; s<pop.numSpecies; s++) numMomentsTot += pop.spec[s].numMoments;

  if ((res == hostOnly) || (res == hostAndDevice))
    mom->ho = alloc_realArray(numMomentsTot*prod_mint(grid.Nx,nDim));  // Allocate on host.
}

void alloc_fourierMoments(const struct fourierGrid grid, const struct population pop, const resource res, struct fourierMoments *momk) {
  mint numMomentsTot = 0;
  for (mint s=0; s<pop.numSpecies; s++) numMomentsTot += pop.spec[s].numMoments;

  if ((res == hostOnly) || (res == hostAndDevice))
    momk->ho = alloc_fourierArray(numMomentsTot*prod_mint(grid.Nekx,nDim));  // Allocate on host.
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
