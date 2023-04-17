/* mugy: field.c
 *
 * Field object and operations to compute it.
 *
 */
#include "mh_field.h"
#include "mh_field_dev.h"
#include "mh_fourier_ho.h"
#include <stdlib.h>  // for malloc.
#include <math.h>  // for pow.

struct mugy_field *mugy_field_alloc() {
  // Allocate the field object.
  struct mugy_field *field = (struct mugy_field *) malloc(sizeof(struct mugy_field));
  return field;
}

void mugy_field_init(struct mugy_field *field, struct mugy_grid *grid, struct mugy_population *pop) {
#ifdef USE_GPU
  enum mugy_resource_mem onResource = MUGY_HOSTDEVICE_MEM;
#else
  enum mugy_resource_mem onResource = MUGY_HOST_MEM;
#endif

  // Allocate Fourier-space potential.
  field->phik = mugy_array_alloc(MUGY_FOURIER, grid->local->fourier->NxTot, onResource);

  // Allocate Fourier-space gyroaveraged potentials, 3 for each species:
  // <phi>, 0.5*hatLap <phi> and (1+0.5*hatLap+hathatLap) <phi>, where
  // <phi> = <J_0>phi is the gyroaveraged potential.
  field->gyrophik = mugy_array_alloc(MUGY_FOURIER, pop->local->numSpecies * 3 * grid->local->fourier->NxTot, onResource);
}

void mugy_field_constop_init(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid, struct mugy_flr *flr) {
  // Pre-compute time independent operators multipying each
  // moment in Poisson equation. Store them in population.

  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_population_species *popL = pop->local,  *popG = pop->global;;

  mint elcIdx = 0;  // Assume electrons are first.
  // Reciprocal of the factor multiplying the potential in the field equation.
  struct mugy_array *rPoiPhikFac = mugy_array_alloc(MUGY_FOURIER, gridL->NxyTot, MUGY_HOST_MEM);
  real delta0e = popG->pars[elcIdx].delta0;
  fourier iimag = 0. + 1.*I;
  for (mint linIdx=0; linIdx<gridL->NxyTot; linIdx++) {
    mint idx[2];  real kx[2];
    mugy_grid_lin2sub_fourier_perp(idx, linIdx, gridL);
    mugy_grid_get_kx_perp(kx, idx, gridL);

    double ky = kx[1];
    fourier *rPoiPhikFac_p = mugy_array_get(rPoiPhikFac, linIdx);

    rPoiPhikFac_p[0] = pow(field->pars.lambdaD,2)*gridL->xperpSq[linIdx]+1.0-delta0e*iimag*ky;
    // Add the (1/tau)*(1-Gamma0Ion) contributions.
    for (mint s=elcIdx+1; s<popG->numSpecies; s++) {  // Assume electrons are first.
      real tau_i = popG->pars[s].tau;
      rPoiPhikFac_p[0] += (1./tau_i)*gridL->xperpSq[linIdx];
    }

    if (MUGY_REAL_MIN < cabs(rPoiPhikFac_p[0]))
      rPoiPhikFac_p[0] = 1./rPoiPhikFac_p[0];
    else
      rPoiPhikFac_p[0] = 0.+0.*I;
  }

  for (mint s=elcIdx+1; s<popL->numSpecies; s++) {
    fourier *denkPoiFac_p  = mugy_population_getMoment_fourier(gridL, popL, s, denIdx, popL->poissonFac->ho);
    fourier *tempkPoiFac_p = mugy_population_getMoment_fourier(gridL, popL, s, tempIdx, popL->poissonFac->ho);

    real tau_i       = popL->pars[s].tau;
    real deltaPerp_i = popL->pars[s].deltaPerp;
    for (mint linIdx=0; linIdx<gridL->NxyTot; linIdx++) {
      fourier *rPoiPhikFac_p = mugy_array_get(rPoiPhikFac, linIdx);

      mint lin1 = s*gridL->NxyTot + linIdx;
      real *avgJ0_p      = mugy_array_get(flr->avgJ0, lin1);
      real *hatLap_p     = mugy_array_get(flr->hatLap, lin1);
      real *hathatLap_p  = mugy_array_get(flr->hathatLap, lin1);

      real hatLapSq      = pow(hatLap_p[0],2);
      real Nb            = 1.+hathatLap_p[0]-0.50*hatLapSq;
      real Db            = 1.+hathatLap_p[0]-0.25*hatLapSq;
      real Sb            = (avgJ0_p[0]*Nb)/Db;
      real hatLapAvgJ0D2 = 0.5*hatLap_p[0]*avgJ0_p[0];

      denkPoiFac_p[0] = rPoiPhikFac_p[0];
      denkPoiFac_p++;

      tempkPoiFac_p[0] = (1./(tau_i*deltaPerp_i*Db))*hatLapAvgJ0D2*rPoiPhikFac_p[0];
      tempkPoiFac_p++;
    };
  };

#ifdef USE_GPU
  mugy_array_copy(popL->poissonFac, popL->poissonFac, MUGY_HOST2DEVICE);
#endif

  mugy_array_free(rPoiPhikFac, MUGY_HOST_MEM);
}

void mugy_field_poisson_solve(struct mugy_field *field, struct mugy_population *pop, struct mugy_grid *grid, mint tstepIdx) {
  // Solve the Poisson equation to obtain phik using the charge density
  // in the time-stepping index 'tstepIdx'
  
#ifdef USE_GPU
  return mugy_field_poisson_solve_dev(field, pop, grid, tstepIdx);
#endif

  struct mugy_population_species *popL = pop->local;
  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_array      *phik  = field->phik;
  struct mugy_array      *momk  = popL->momk[tstepIdx];

  mugy_array_zero(phik, MUGY_HOST_MEM);

  for (mint s=0; s<popL->numSpecies; s++) {
    for (mint m=0; m<popL->pars[s].numMoments; m++) {
      fourier *poissonFac_p = mugy_population_getMoment_fourier(gridL, popL, s, m, popL->poissonFac->ho);
      fourier *momk_p       = mugy_population_getMoment_fourier(gridL, popL, s, m, momk->ho);
      for (mint linIdx=0; linIdx<gridL->NxTot; linIdx++) {
        fourier *phik_p = mugy_array_get(phik, linIdx);

	phik_p[0] += poissonFac_p[0]*momk_p[0];

	poissonFac_p++;  momk_p++;
      }
    }
  }
}

void mugy_field_free(struct mugy_field *field) {
#ifdef USE_GPU
  enum mugy_resource_mem onResource = MUGY_HOSTDEVICE_MEM;
#else
  enum mugy_resource_mem onResource = MUGY_HOST_MEM;
#endif

  // Free potentials.
  mugy_array_free(field->phik, onResource); 
  mugy_array_free(field->gyrophik, onResource); 

  // Free the field object.
  free(field);
}
