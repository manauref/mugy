/* mugy: linear.c
 *
 * Pre-compute the linear operators.
 *
 */
#include "mh_linear.h"
#include "mh_fourier_ho.h"

void mugy_linear_constop_init(struct mugy_population *pop, struct mugy_grid *grid, struct mugy_flr *flr) {
  // Pre-compute the time independent linear operators in each equation.

  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_population_species *popL = pop->local;

  fourier iimag = 0. + 1.*I;

  // Linear operators acting on the potential.
  for (mint s=0; s<popL->numSpecies; s++) {
    real charge = popL->pars[s].qCharge;  real omd       = popL->pars[s].omd;
    real tau    = popL->pars[s].tau;      real omSt      = popL->pars[s].omSt;
    real eta    = popL->pars[s].eta;      real deltaPerp = popL->pars[s].deltaPerp;

    fourier *denkOp_p  = mugy_population_getMoment_fourier(gridL, popL, s, denIdx, popL->linOpPhi->ho);
    fourier *tempkOp_p = mugy_population_getMoment_fourier(gridL, popL, s, tempIdx, popL->linOpPhi->ho);

    for (mint linIdx=0; linIdx<gridL->NxTot; linIdx++) {
      mint idx[nDim];  real kx[nDim];
      mugy_grid_lin2sub_fourier(idx, linIdx, gridL);
      mugy_grid_get_kx(kx, idx, gridL);

      mint lin1 = s*gridL->NxyTot + linIdx;
      real *hatLap_p    = mugy_array_get(flr->hatLap   , lin1);
      real *hathatLap_p = mugy_array_get(flr->hathatLap, lin1);

      double ky = kx[1];
      fourier iOmd  = -charge*tau*iimag*omd*ky;
      fourier iOmSt =        -tau*iimag*omSt*ky;
      
      denkOp_p[0] = (1./tau)*((1.+eta*0.5*hatLap_p[0])*iOmSt-(2.+0.5*hatLap_p[0])*iOmd);

      tempkOp_p[0] = deltaPerp*(((1.+eta)*(1.+0.5*hatLap_p[0])
                    +eta*hathatLap_p[0])*iOmSt-charge*(3.+1.5*hatLap_p[0]+hathatLap_p[0])*iOmd);
    }
  }

  // Linear operators acting on the moments. We save these as a column-major order
  // matrix for each species in case we wish to apply them using cuBLAS.
  mint sOff = 0;
  for (mint s=0; s<popL->numSpecies; s++) {
    real charge = popL->pars[s].qCharge;  real omd   = popL->pars[s].omd;
    real tau    = popL->pars[s].tau;      real omSt  = popL->pars[s].omSt;
    real *alpha = popL->pars[s].alpha;    real delta = popL->pars[s].delta;
    real *nu    = popL->pars[s].nu;
  
    for (mint linIdx=0; linIdx<gridL->NxTot; linIdx++) {
      mint idx[nDim];  real kx[nDim];
      mugy_grid_lin2sub_fourier(idx, linIdx, gridL);
      mugy_grid_get_kx(kx, idx, gridL);

      // Convert the 3D linIdx to a 2D (perp) linIdx to index into xperpSq.
      mint idxperp[2];
      for (mint d=0; d<2; d++) idxperp[d] = idx[d];
      mint linIdxperp = mugy_grid_sub2lin_perp_fourier(idxperp, gridL);
      real kperpSq = gridL->xperpSq[linIdxperp];

      double ky = kx[1];
      fourier iOmd  = -charge*tau*iimag*omd*ky;
      fourier iOmSt =        -tau*iimag*omSt*ky;
      
      // IMPORTANT: Notice the column-major order of these operators.
      fourier *den_den_p   = mugy_array_get(popL->linOpMom, linIdx+sOff+0*gridL->NxTot);
      fourier *temp_den_p  = mugy_array_get(popL->linOpMom, linIdx+sOff+1*gridL->NxTot);
      fourier *den_temp_p  = mugy_array_get(popL->linOpMom, linIdx+sOff+2*gridL->NxTot);
      fourier *temp_temp_p = mugy_array_get(popL->linOpMom, linIdx+sOff+3*gridL->NxTot);

      den_den_p[0]  = -(1./tau)*iOmd*tau*delta-alpha[0]-nu[0]*kperpSq;
      den_temp_p[0] = -(1./tau)*iOmd;

      temp_den_p[0]  = 0.;
      temp_temp_p[0] = -alpha[1]-nu[1]*kperpSq;
    }
    sOff += popL->pars[s].numMoments * popL->pars[s].numMoments * gridL->NxTot;
  }

#ifdef USE_GPU
  mugy_array_copy(popL->linOpMom, popL->linOpMom, MUGY_HOST2DEVICE);
  mugy_array_copy(popL->linOpPhi, popL->linOpPhi, MUGY_HOST2DEVICE);
#endif
}
