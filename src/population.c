/* mugy: population.c
 *
 * Methods to work with data/information of each species in the population object.
 *
 */

#include "mh_population.h"

real* getMoment_real(struct mugy_realGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, real *momIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in mom.
  real* ptrOut = momIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NxTot;
}

fourier* getMoment_fourier(struct mugy_fourierGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, fourier *momkIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
  fourier* ptrOut = momkIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NekxTot;
}
