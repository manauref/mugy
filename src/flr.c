/* mugy: flr.c
 *
 * Create various FLR (finite Larmor radius) operators.
 *
 */
#include "mh_flr.h"
#include "mh_fourier_ho.h"
#include "mh_alloc.h"
#include "mh_array.h"
#include "mh_bessel.h"
#include <math.h>
#include <float.h> // For FLT_MIN, DBL_MIN.

#ifdef USE_SINGLE_PRECISION
const real REAL_MIN = FLT_MIN;
#else
const real REAL_MIN = DBL_MIN;
#endif

void mugy_flr_init(struct mugy_population *pop, struct mugy_grid *grids, struct mugy_field *field, struct mugy_ioManager *ioman) {
  // Initialize the FLR factors, store them in population.

  struct mugy_fourierGrid *grid = &grids->local.deal;
  struct mugy_pop *popL = &pop->local;

  // Compute the argument of the Bessel functions, b=(kperp*rho)^2.
  struct mugy_array *bFLR = mugy_array_alloc(real_enum, popL->numSpecies * grid->NekxyTot, hostMem);
  for (mint linIdx=0; linIdx<grid->NekxyTot; linIdx++) {
    real kperp = sqrt(grid->kperpSq[linIdx]);
    for (mint s=0; s<popL->numSpecies; s++) {
      real muMass   = pop->global.pars[s].muMass;
      real tau      = pop->global.pars[s].tau;
      real kperprho = sqrt(tau)*kperp/muMass;

      real *bFLR_p = mugy_array_get(bFLR, s*grid->NekxyTot + linIdx);
      bFLR_p[0]    = pow(kperprho,2);
    }
  }

#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Allocate space for the FLR operators inside Poisson brackets, 3 for
  // each species: <J_0>=Gamma_0^{1/2}, 0.5*hatLap <J_0>, (1+0.5*hatLap+hathatLap) <J_0>.
  popL->pbFLRop = mugy_array_alloc(real_enum, popL->numSpecies * 3 * grid->NekxyTot, onResource);

  struct mugy_array *Gamma0    = mugy_array_alloc(real_enum, popL->numSpecies * grid->NekxyTot, hostMem);
  struct mugy_array *hatLap    = mugy_array_alloc(real_enum, popL->numSpecies * grid->NekxyTot, hostMem);
  struct mugy_array *hathatLap = mugy_array_alloc(real_enum, popL->numSpecies * grid->NekxyTot, hostMem);
  struct mugy_array *poissonDb = mugy_array_alloc(real_enum, popL->numSpecies * grid->NekxyTot, hostMem);
  struct mugy_array *poissonSb = mugy_array_alloc(real_enum, popL->numSpecies * grid->NekxyTot, hostMem);
  real *pbFLRop_p = popL->pbFLRop->ho;
  for (mint s=0; s<popL->numSpecies; s++) {
    for (mint linIdx=0; linIdx<grid->NekxyTot; linIdx++) {
      mint lin1 = s*grid->NekxyTot + linIdx;

      real *bFLR_p      = mugy_array_get(bFLR     , lin1);
      real *Gamma0_p    = mugy_array_get(Gamma0   , lin1);
      real *hatLap_p    = mugy_array_get(hatLap   , lin1);
      real *hathatLap_p = mugy_array_get(hathatLap, lin1);
      real *poissonDb_p = mugy_array_get(poissonDb, lin1);
      real *poissonSb_p = mugy_array_get(poissonSb, lin1);

      real b = bFLR_p[0];

      Gamma0_p[0]    = mugy_bessel_I0_scaled(b);
      real GammaRat  = mugy_bessel_I1_scaled(b)/Gamma0_p[0];
      hatLap_p[0]    = b*(GammaRat-1.);
      hathatLap_p[0] = b*((0.5*GammaRat-1.)-0.25*b*(3.+GammaRat)*(GammaRat-1.));

      mint lin3 = s*3*grid->NekxyTot + linIdx;
      real avgJ0 = sqrt(Gamma0_p[0]);  // <J_0>.
      pbFLRop_p[lin3+0*grid->NekxyTot] = avgJ0; 
      pbFLRop_p[lin3+1*grid->NekxyTot] = 0.5*hatLap_p[0]*avgJ0;
      pbFLRop_p[lin3+2*grid->NekxyTot] = (1.+0.5*hatLap_p[0]+hathatLap_p[0])*avgJ0;

      real hatLapSq  = pow(hatLap_p[0],2);
      real Nb        = 1.+hathatLap_p[0]-0.50*hatLapSq;
      poissonDb_p[0] = 1.+hathatLap_p[0]-0.25*hatLapSq;
      poissonSb_p[0] = (avgJ0*Nb)/poissonDb_p[0];
    }
  }

#ifdef USE_GPU
  mugy_array_copy(popL->pbFLRop, popL->pbFLRop, host2device);
#endif

//    struct mugy_ad_file *fhr = mugy_io_create_population_perp_file(ioman, "pbFLRop", *grids, *pop, real_enum, fourier_enum, 3, 0);
//    mugy_io_write_mugy_array(NULL, "pbFLRop", fhr, popL->pbFLRop);
//    mugy_io_close_file(fhr);
//
//    struct mugy_ad_file *fhs = mugy_io_create_population_perp_file(ioman, "poissonDb", *grids, *pop, real_enum, fourier_enum, 1, 0);
//    mugy_io_write_mugy_array(NULL, "poissonDb", fhs, poissonDb);
//    mugy_io_close_file(fhs);
//    struct mugy_ad_file *fht = mugy_io_create_population_perp_file(ioman, "poissonSb", *grids, *pop, real_enum, fourier_enum, 1, 0);
//    mugy_io_write_mugy_array(NULL, "poissonSb", fht, poissonSb);
//    mugy_io_close_file(fht);


  mint elcIdx = 0;  // Assume electrons are first.
  // Reciprocal of the factor multiplying the potential in the field equation.
  struct mugy_array *rPoiPhikFac = mugy_array_alloc(fourier_enum, grid->NekxyTot, hostMem);
  real delta0e = pop->global.pars[elcIdx].delta0;
  fourier iimag = 0. + 1.*I;
  for (mint linIdx=0; linIdx<grid->NekxyTot; linIdx++) {
    mint idx[2];  real kx[2];
    mugy_grid_lin2sub_fourier_perp(idx, linIdx, *grid);
    mugy_grid_get_kx_perp(kx, idx, *grid);

    double ky = kx[1];
    fourier *rPoiPhikFac_p = mugy_array_get(rPoiPhikFac, linIdx);

    rPoiPhikFac_p[0] = pow(field->pars.lambdaD,2)*grid->kperpSq[linIdx]+1.0-delta0e*iimag*ky;
    // Add the (1/tau)*(1-Gamma0Ion) contributions.
    for (mint s=elcIdx+1; s<pop->global.numSpecies; s++) {  // Assume electrons are first.
      real tau_i = pop->global.pars[s].tau;
      rPoiPhikFac_p[0] += (1./tau_i)*grid->kperpSq[linIdx]; 
    }

    if (REAL_MIN < cabs(rPoiPhikFac_p[0]))
      rPoiPhikFac_p[0] = 1./rPoiPhikFac_p[0]; 
    else
      rPoiPhikFac_p[0] = 0.+0.*I;
  }

  // Pre-compute the factors multipying each moment in the Poisson equation.
  popL->poissonFac = mugy_array_alloc(fourier_enum, popL->numMomentsTot*grid->NekxyTot, onResource);

  for (mint s=elcIdx+1; s<popL->numSpecies; s++) {
    fourier *denkPoiFac_p  = mugy_population_getMoment_fourier(*grid, *popL, s, denIdx, popL->poissonFac->ho);
    fourier *tempkPoiFac_p = mugy_population_getMoment_fourier(*grid, *popL, s, tempIdx, popL->poissonFac->ho);

    real tau_i       = popL->pars[s].tau;
    real deltaPerp_i = popL->pars[s].deltaPerp;
    for (mint linIdx=0; linIdx<grid->NekxyTot; linIdx++) {
      fourier *rPoiPhikFac_p = mugy_array_get(rPoiPhikFac, linIdx);

      mint lin1 = s*grid->NekxyTot + linIdx;
      real *poissonDb_p = mugy_array_get(poissonDb, lin1);

      mint lin3 = s*3*grid->NekxyTot + linIdx;
      real *pbFLRop_p = mugy_array_get(popL->pbFLRop, lin3+1*grid->NekxyTot);; 
      real hatLapAvgJ0D2 = pbFLRop_p[0];

      denkPoiFac_p[0] = rPoiPhikFac_p[0];
      denkPoiFac_p++;

      tempkPoiFac_p[0] = (1./(tau_i*deltaPerp_i*poissonDb_p[0]))*hatLapAvgJ0D2*rPoiPhikFac_p[0];
      tempkPoiFac_p++;
    };
  };

#ifdef USE_GPU
  mugy_array_copy(popL->poissonFac, popL->poissonFac, host2device);
#endif

  mugy_array_free(rPoiPhikFac, hostMem);
  mugy_array_free(poissonSb, hostMem);
  mugy_array_free(poissonDb, hostMem);
  mugy_array_free(hathatLap, hostMem);
  mugy_array_free(hatLap, hostMem);
  mugy_array_free(Gamma0, hostMem);
  mugy_array_free(bFLR, hostMem);
}

