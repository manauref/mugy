/* mugy: nonlinear.c
 *
 * Infrastructure to compute the nonlinear terms (i.e. Poisson brackets).
 *
 */
#include "mh_nonlinear.h"

void mugy_nonlinear_constop_init(struct mugy_population *pop, struct mugy_grid *grid, struct mugy_flr *flr) {
  // Pre-compute the time independent operators that multiply the potential phi inside the Poisson brackets:
  // <J_0>=Gamma_0^{1/2}, 0.5*hatLap <J_0>, (1+0.5*hatLap+hathatLap) <J_0>.

  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_population_species *popL = pop->local;

  real *pbFLRop_p = popL->pbFLRop->ho;
  for (mint s=0; s<popL->numSpecies; s++) {
    for (mint linIdx=0; linIdx<gridL->NxyTot; linIdx++) {
      mint lin1 = s*gridL->NxyTot + linIdx;
      real *avgJ0_p     = mugy_array_get(flr->avgJ0    , lin1);
      real *hatLap_p    = mugy_array_get(flr->hatLap   , lin1);
      real *hathatLap_p = mugy_array_get(flr->hathatLap, lin1);

      mint lin3 = s*3*gridL->NxyTot + linIdx;
      pbFLRop_p[lin3+0*gridL->NxyTot] = avgJ0_p[0];
      pbFLRop_p[lin3+1*gridL->NxyTot] = 0.5*hatLap_p[0]*avgJ0_p[0];
      pbFLRop_p[lin3+2*gridL->NxyTot] = (1.+0.5*hatLap_p[0]+hathatLap_p[0])*avgJ0_p[0];
    }
  }

#ifdef USE_GPU
  mugy_array_copy(popL->pbFLRop, popL->pbFLRop, MUGY_HOST2DEVICE);
#endif
}
