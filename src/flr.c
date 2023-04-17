/* mugy: flr.c
 *
 * Create various FLR (finite Larmor radius) operators.
 *
 */
#include "mh_flr.h"
#include "mh_bessel.h"
#include <math.h>
#include <stdlib.h>  // for malloc.

struct mugy_flr *mugy_flr_init(struct mugy_population *pop, struct mugy_grid *grid) {
  // Initialize the FLR operators on the host.
  struct mugy_flr* flr = (struct mugy_flr*) malloc(sizeof(struct mugy_flr)); 

  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_population_species *popL = pop->local;

  // Compute the argument of the Bessel functions, b=(kperp*rho)^2.
  struct mugy_array *bFLR = mugy_array_alloc(MUGY_REAL, popL->numSpecies * gridL->NxyTot, MUGY_HOST_MEM);
  for (mint linIdx=0; linIdx<gridL->NxyTot; linIdx++) {
    real kperp = sqrt(gridL->xperpSq[linIdx]);

    for (mint s=0; s<popL->numSpecies; s++) {
      real muMass   = pop->global->pars[s].muMass;
      real tau      = pop->global->pars[s].tau;
      real kperprho = sqrt(tau)*kperp/muMass;

      real *bFLR_p = mugy_array_get(bFLR, s*gridL->NxyTot + linIdx);
      bFLR_p[0]    = pow(kperprho,2);
    }
  }

  // Compute the FLR operators that other (FLR) operators are built upon:
  //   Gamma0: exponentially weighted 0th modified Bessel function of the 1st kind.
  //   <J_0>=Gamma0^{1/2}.
  //   \hat{nabla}_perp.
  //   \hat{hat{nabla}}_perp.
  // These only depend on kx-ky, so we defined them on the perp plane.
  flr->Gamma0    = mugy_array_alloc(MUGY_REAL, popL->numSpecies * gridL->NxyTot, MUGY_HOST_MEM);
  flr->avgJ0     = mugy_array_alloc(MUGY_REAL, popL->numSpecies * gridL->NxyTot, MUGY_HOST_MEM);
  flr->hatLap    = mugy_array_alloc(MUGY_REAL, popL->numSpecies * gridL->NxyTot, MUGY_HOST_MEM);
  flr->hathatLap = mugy_array_alloc(MUGY_REAL, popL->numSpecies * gridL->NxyTot, MUGY_HOST_MEM);
  for (mint linIdx=0; linIdx<popL->numSpecies*gridL->NxyTot; linIdx++) {
    real *bFLR_p      = mugy_array_get(bFLR, linIdx);
    real *Gamma0_p    = mugy_array_get(flr->Gamma0   , linIdx);
    real *avgJ0_p     = mugy_array_get(flr->avgJ0    , linIdx);
    real *hatLap_p    = mugy_array_get(flr->hatLap   , linIdx);
    real *hathatLap_p = mugy_array_get(flr->hathatLap, linIdx);

    real b = bFLR_p[0];

    Gamma0_p[0]    = mugy_bessel_I0_scaled(b);
    avgJ0_p[0]     = sqrt(Gamma0_p[0]);
    real GammaRat  = mugy_bessel_I1_scaled(b)/Gamma0_p[0];
    hatLap_p[0]    = b*(GammaRat-1.);
    hathatLap_p[0] = b*((0.5*GammaRat-1.)-0.25*b*(3.+GammaRat)*(GammaRat-1.));
  }

  mugy_array_free(bFLR, MUGY_HOST_MEM);

  return flr;
}

void mugy_flr_free(struct mugy_flr *flr) {
  mugy_array_free(flr->Gamma0   , MUGY_HOST_MEM);
  mugy_array_free(flr->avgJ0    , MUGY_HOST_MEM);
  mugy_array_free(flr->hatLap   , MUGY_HOST_MEM);
  mugy_array_free(flr->hathatLap, MUGY_HOST_MEM);
  free(flr);
}
