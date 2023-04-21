/* mugy: dydt.c
 *
 * Compute the time rate of change of the moment vector.
 *
 */
#include "mh_dydt.h"
#include "mh_dydt_dev.h"
#include "mh_fourier_ho.h"

void mugy_dydt_phikFLR(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid) {
  // Compute <J_0>phik, 0.5*hatLap<J_0>phik, 1+0.5*hatLap+hathatLap)<J_0>phik.

#if USE_GPU
//  return mugy_dydt_phikFLR(pop, field, grid);
#endif

  struct mugy_population_species *popL = pop->local;
  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_array *phik = field->phik;

  for (mint s=0; s<popL->numSpecies; s++) {
    for (unsigned long linIdx=0; linIdx<gridL->NxTot; linIdx++) {
      fourier *phik_p = mugy_array_get(phik, linIdx);

      unsigned long lin3 = s*3*gridL->NxTot + linIdx;
      fourier *phikFLR0_p = mugy_array_get(popL->phikFLR, lin3+0*gridL->NxTot);
      fourier *phikFLR1_p = mugy_array_get(popL->phikFLR, lin3+1*gridL->NxTot);
      fourier *phikFLR2_p = mugy_array_get(popL->phikFLR, lin3+2*gridL->NxTot);

      // Convert the 3D linIdx to a 2D (perp) linIdx to index FLR ops.
      mint idx[nDim], idxperp[2];
      mugy_grid_lin2sub_fourier(idx, linIdx, gridL);
      for (mint d=0; d<2; d++) idxperp[d] = idx[d];
      unsigned long linIdxperp = mugy_grid_sub2lin_perp_fourier(idxperp, gridL);

      unsigned long lin3perp = s*3*gridL->NxyTot + linIdxperp;
      real *flrOp0_p = mugy_array_get(popL->pbFLRop, lin3perp+0*gridL->NxyTot);
      real *flrOp1_p = mugy_array_get(popL->pbFLRop, lin3perp+1*gridL->NxyTot);
      real *flrOp2_p = mugy_array_get(popL->pbFLRop, lin3perp+2*gridL->NxyTot);

      phikFLR0_p[0] = flrOp0_p[0]*phik_p[0];
      phikFLR1_p[0] = flrOp1_p[0]*phik_p[0];
      phikFLR2_p[0] = flrOp2_p[0]*phik_p[0];
    }
  }
}

void mugy_dydt_linear(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid) {
  // Compute the time rate of change dy/dt, where y is the vector of moments, due to linear terms.

//#if USE_GPU
//  return mugy_dydt_linear_dev(outIdx, inIdx, time, pop, field, grid);
//#endif

  struct mugy_population_species *popL = pop->local;
  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_array *momkDot = popL->momk[outIdx];
  struct mugy_array *momk    = popL->momk[inIdx];

  // Add the linear terms.
  mint NxTot = gridL->NxTot;
  mint sOff = 0;
  for (mint s=0; s<popL->numSpecies; s++) {
    mint numMom = popL->pars[s].numMoments;
    for (mint m=0; m<numMom; m++) {
      // Add contributions to the time rate of change of the m-th moment.
      for (unsigned long linIdx=0; linIdx<NxTot; linIdx++) {
        fourier *momkDot_p = mugy_population_getMoment_fourier(gridL, popL, s, m, momkDot->ho);

        // Linear operator acting on the n-th moment. 
        for (mint n=0; n<m; n++) {
          fourier *momk_p     = mugy_population_getMoment_fourier(gridL, popL, s, n, momk->ho);
          fourier *linOpMom_p = mugy_array_get(popL->linOpMom, sOff+(m*numMom+n)*NxTot+linIdx);
          momkDot_p[0] += linOpMom_p[0]*momk_p[0];
        }
        
        // Linear operator acting on the potential. 
        fourier *linOpPhi_p = mugy_population_getMoment_fourier(gridL, popL, s, m, popL->linOpPhi->ho);
        unsigned long lin3 = s*3*gridL->NxTot + linIdx;
        fourier *phikGyroAvg_p = mugy_array_get(popL->phikFLR, lin3);
        momkDot_p[0] += linOpPhi_p[0]*phikGyroAvg_p[0];
      }
    }
    sOff += numMom*numMom*NxTot;
  }

}

void mugy_dydt(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid, struct mugy_fft *fft) {
  // Compute the time rate of change dy/dt, where y is the vector of moments.

#if USE_GPU
  return mugy_dydt_dev(outIdx, inIdx, time, pop, field, grid, fft);
#endif

  struct mugy_population_species *popL = pop->local;
  struct mugy_array *momkDot = popL->momk[outIdx];
  mugy_array_zero(momkDot, MUGY_HOST_MEM);

  // Compute the gyroaveraged potentials.
  mugy_dydt_phikFLR(pop, field, grid); 

  // Apply linear terms.
  mugy_dydt_linear(outIdx, inIdx, time, pop, field, grid); 
}
