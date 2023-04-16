/* mugy: field.c
 *
 * Field object and operations to compute it.
 *
 */
#include "mh_field.h"
#include "mh_field_dev.h"
#include "mh_fourier_ho.h"
#include <stdlib.h>  // for malloc.

struct mugy_field *mugy_field_alloc() {
  // Allocate the field object.
  struct mugy_field *field = (struct mugy_field *) malloc(sizeof(struct mugy_field));
  return field;
}

void mugy_field_init(struct mugy_field *field, struct mugy_grid *grid, struct mugy_population *pop) {
#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Allocate Fourier-space potential.
  field->phik = mugy_array_alloc(fourier_enum, grid->local.deal.NekxTot, onResource);

  // Allocate Fourier-space gyroaveraged potentials, 3 for each species:
  // <phi>, 0.5*hatLap <phi> and (1+0.5*hatLap+hathatLap) <phi>, where
  // <phi> = <J_0>phi is the gyroaveraged potential.
  field->gyrophik = mugy_array_alloc(fourier_enum, pop->local.numSpecies * 3 * grid->local.deal.NekxTot, onResource);
}

void mugy_field_free(struct mugy_field *field) {
#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Free potentials.
  mugy_array_free(field->phik, onResource); 
  mugy_array_free(field->gyrophik, onResource); 

  // Free the field object.
  free(field);
}

void mugy_field_poisson_solve(struct mugy_field *field, struct mugy_population *pop, struct mugy_grid *grid, mint tstepIdx) {
  // Solve the Poisson equation to obtain phik using the charge density
  // in the time-stepping index 'tstepIdx'
  
#ifdef USE_GPU
  return mugy_field_poisson_solve_dev(field, pop, grid, tstepIdx);
#endif

  struct mugy_fourierGrid *gridL = &grid->local.deal;
  struct mugy_pop         *popL  = &pop->local;
  struct mugy_array       *phik  = field->phik;
  struct mugy_array       *momk  = popL->momk[tstepIdx];

  mugy_array_zero(phik, hostMem);

  for (mint s=0; s<popL->numSpecies; s++) {
    for (mint m=0; m<popL->pars[s].numMoments; m++) {
      fourier *poissonFac_p = mugy_population_getMoment_fourier(*gridL, *popL, s, m, popL->poissonFac->ho);
      fourier *momk_p       = mugy_population_getMoment_fourier(*gridL, *popL, s, m, momk->ho);
      for (mint linIdx=0; linIdx<gridL->NekxTot; linIdx++) {
        fourier *phik_p = mugy_array_get(phik, linIdx);

	phik_p[0] += poissonFac_p[0]*momk_p[0];

	poissonFac_p++;  momk_p++;
      }
    }
  }
}
