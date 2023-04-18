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

  // Allocate space for local and global populations.
  pop->local  = (struct mugy_population_species *) malloc(sizeof(struct mugy_population_species));
  pop->global = (struct mugy_population_species *) malloc(sizeof(struct mugy_population_species));

  // Set pointers to large global arrays to NULL in global population (don't store them).
  pop->global->momk = NULL;
  pop->global->pbFLRop = NULL;
  pop->global->poissonFac = NULL;

  return pop;
}

// Functions that allocate moment vectors.
struct mugy_array *mugy_population_alloc_realMoments(struct mugy_grid_basic *grid,
  struct mugy_population_species *pop, enum mugy_resource_mem res) {

  mint nelem = pop->numMomentsTot*grid->NxTot;
  struct mugy_array *mom = mugy_array_alloc(MUGY_REAL, nelem, res);
  return mom;
}
struct mugy_array *mugy_population_alloc_fourierMoments(struct mugy_grid_basic *grid,
  struct mugy_population_species *pop, enum mugy_resource_mem res) {

  mint nelem = pop->numMomentsTot*grid->NxTot;
  struct mugy_array *momk = mugy_array_alloc(MUGY_FOURIER, nelem, res);
  return momk;
}

void mugy_population_alloc_local(struct mugy_population_species *popL, struct mugy_grid_chart *gridL) {
  // Allocate various moments needed.
#ifdef USE_GPU
  enum mugy_resource_mem onResource = MUGY_HOSTDEVICE_MEM;
#else
  enum mugy_resource_mem onResource = MUGY_HOST_MEM;
#endif

  // Allocate moments vector needed for time stepping.
  popL->momk = (struct mugy_array**) calloc(TIME_STEPPER_NUM_FIELDS, sizeof(struct mugy_array *));
  for (mint s=0; s<TIME_STEPPER_NUM_FIELDS; s++)
    popL->momk[s] = mugy_population_alloc_fourierMoments(gridL->fourier, popL, onResource);

  // Allocate space for linear operators.
  popL->linOpPhi = mugy_array_alloc(MUGY_FOURIER, popL->numMomentsTot*gridL->fourier->NxTot, onResource);
  mint numOps = 0;
  for (mint s=0; s<popL->numSpecies; s++) numOps += popL->pars[s].numMoments * popL->pars[s].numMoments;
  popL->linOpMom = mugy_array_alloc(MUGY_FOURIER, numOps*gridL->fourier->NxTot, onResource);

  // Allocate space for factors multipying each moment in the Poisson equation.
  popL->poissonFac = mugy_array_alloc(MUGY_FOURIER, popL->numMomentsTot*gridL->fourier->NxyTot, onResource);

  // Allocate space for the FLR operators inside Poisson brackets, 3 for
  // each species: <J_0>=Gamma_0^{1/2}, 0.5*hatLap <J_0>, (1+0.5*hatLap+hathatLap) <J_0>.
  popL->pbFLRop = mugy_array_alloc(MUGY_REAL, popL->numSpecies * 3 * gridL->fourier->NxyTot, onResource);
}

real* mugy_population_getMoment_real(struct mugy_grid_basic *grid, struct mugy_population_species *pop, mint sIdx, mint momIdx, real *momIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in mom.
  real* ptrOut = momIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop->pars[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid->NxTot;
}

void* mugy_population_getMoment_fourier(struct mugy_grid_basic *grid, struct mugy_population_species *pop, mint sIdx, mint momIdx, void *momkIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
  fourier* ptrOut = (fourier *)momkIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop->pars[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid->NxTot;
}

void mugy_population_free(struct mugy_population *pop) {
  // Deallocate memory used in species struct.
  for (mint i=0; i<2; i++) {
    struct mugy_population_species *spec = i==0 ? pop->global : pop->local;
    for (mint s=0; s<spec->numSpecies; s++) {
      free(spec->pars[s].alpha);
      free(spec->pars[s].nu);
      free(spec->pars[s].hDiffOrder);
      free(spec->pars[s].hDiff);
      free(spec->pars[s].kDiffMin);
      free(spec->pars[s].initAux);
    }
    free(spec->pars);
  }

#ifdef USE_GPU
  enum mugy_resource_mem onResource = MUGY_HOSTDEVICE_MEM;
#else
  enum mugy_resource_mem onResource = MUGY_HOST_MEM;
#endif

  // Free moments vector.
  for (mint s=0; s<TIME_STEPPER_NUM_FIELDS; s++)
    mugy_array_free(pop->local->momk[s], onResource);
  free(pop->local->momk);

  // Free linear operators.
  mugy_array_free(pop->local->linOpPhi, onResource);
  mugy_array_free(pop->local->linOpMom, onResource);

  // Free factors multiplying moments in Poisson equation.
  mugy_array_free(pop->local->poissonFac, onResource);

  // Free Poisson bracket FLR operators.
  mugy_array_free(pop->local->pbFLRop, onResource);

  free(pop->local);
  free(pop->global);
  free(pop);
}
