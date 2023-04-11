/* mugy: population.c
 *
 * Methods to work with data/information of each species in the population object.
 *
 */

#include "mh_population.h"
#include "mh_fourier_ho.h"

// Functions that allocate moment vectors.
void alloc_realMoments(struct mugy_array *mom, const struct mugy_realGrid grid,
  const struct mugy_population pop, enum resource_mem res) {
  mint nelem = pop.numMomentsTot*grid.NxTot;
  mugy_array_alloc(mom, real_enum, nelem, res);
}
void alloc_fourierMoments(struct mugy_array *momk, const struct mugy_fourierGrid grid,
  const struct mugy_population pop, enum resource_mem res) {
  mint nelem = pop.numMomentsTot*grid.NekxTot;
  mugy_array_alloc(momk, fourier_enum, nelem, res);
}

real* getMoment_real(struct mugy_realGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, real *momIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in mom.
  real* ptrOut = momIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NxTot;
}

void* getMoment_fourier(struct mugy_fourierGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, void *momkIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
  fourier* ptrOut = (fourier *)momkIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NekxTot;
}
