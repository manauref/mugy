/* mugy: population.c
 *
 * Methods to work with data/information of each species in the population object.
 *
 */

#include "mh_population.h"
#include "mh_fourier_ho.h"
#include <stdlib.h>  // for malloc.

// Allocate mugy_grid object.
struct mugy_population *mugy_population_alloc() {
  struct mugy_population *pop = (struct mugy_population *) malloc(sizeof(struct mugy_population));

  // Set the moments pointer in the global population to NULL (we don't store global moments).
  pop->global.momk = NULL;

  return pop;
}

// Functions that allocate moment vectors.
struct mugy_array *mugy_population_alloc_realMoments(const struct mugy_realGrid grid,
  const struct mugy_pop pop, enum mugy_resource_mem res) {

  mint nelem = pop.numMomentsTot*grid.NxTot;
  struct mugy_array *mom = mugy_array_alloc(MUGY_REAL, nelem, res);
  return mom;
}
struct mugy_array *mugy_population_alloc_fourierMoments(const struct mugy_fourierGrid grid,
  const struct mugy_pop pop, enum mugy_resource_mem res) {

  mint nelem = pop.numMomentsTot*grid.NekxTot;
  struct mugy_array *momk = mugy_array_alloc(MUGY_FOURIER, nelem, res);
  return momk;
}

void mugy_population_alloc_moments(struct mugy_population *pop, struct mugy_grid *grid) {
  // Allocate various moments needed.
#ifdef USE_GPU
  enum mugy_resource_mem onResource = MUGY_HOSTDEVICE_MEM;
#else
  enum mugy_resource_mem onResource = MUGY_HOST_MEM;
#endif

  // Allocate moments vector needed for time stepping.
  pop->local.momk = (struct mugy_array**) calloc(TIME_STEPPER_NUM_FIELDS, sizeof(struct mugy_array *));
  for (mint s=0; s<TIME_STEPPER_NUM_FIELDS; s++)
    pop->local.momk[s] = mugy_population_alloc_fourierMoments(grid->local.deal, pop->local, onResource);

}

real* mugy_population_getMoment_real(struct mugy_realGrid grid, struct mugy_pop pop, mint sIdx, mint momIdx, real *momIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in mom.
  real* ptrOut = momIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.pars[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NxTot;
}

void* mugy_population_getMoment_fourier(struct mugy_fourierGrid grid, struct mugy_pop pop, mint sIdx, mint momIdx, void *momkIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
  fourier* ptrOut = (fourier *)momkIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.pars[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NekxTot;
}

void mugy_population_free(struct mugy_population *pop) {
  // Deallocate memory used in species struct.
  for (mint i=0; i<2; i++) {
    struct mugy_pop *popp = i==0 ? &pop->global : &pop->local;
    for (mint s=0; s<popp->numSpecies; s++) {
      free(popp->pars[s].alpha);
      free(popp->pars[s].nu);
      free(popp->pars[s].hDiffOrder);
      free(popp->pars[s].hDiff);
      free(popp->pars[s].kDiffMin);
      free(popp->pars[s].initAux);
    }
    free(popp->pars);
  }

#ifdef USE_GPU
  enum mugy_resource_mem onResource = MUGY_HOSTDEVICE_MEM;
#else
  enum mugy_resource_mem onResource = MUGY_HOST_MEM;
#endif

  // Free moments vector.
  for (mint s=0; s<TIME_STEPPER_NUM_FIELDS; s++)
    mugy_array_free(pop->local.momk[s], onResource);
  free(pop->local.momk);

  // Free Poisson bracket FLR operators.
  mugy_array_free(pop->local.pbFLRop, onResource);

  // Free factors multiplying moments in Poisson equation.
  mugy_array_free(pop->local.poissonFac, onResource);
}
