/* mugy: mh_population.h
 *
 * Population object holding the data and information about
 * each species (including its moments).
 *
 */
#ifndef MUGY_POPULATION
#define MUGY_POPULATION

#include "mh_data.h"
#include "mh_grid.h"

struct mugy_species {
  mint numMoments;   // Number of moments.
  real qCharge;      // Charge.
  real muMass;       // sqrt of the mass.
  real tau;          // Temperature. 
  real omSt;         // omega_star.
  real omd;          // omega_d.
  real delta;        // (Tpar + Tperp)/T.
  real deltaPerp;    // Tperp/T.
  real eta;          // L_n/L_{Tperp}.
  real *alpha;       // Damping coefficient.
  real *nu;          // Diffusion coefficient.
  real delta0;       // i-delta non-adiabaticity parameter.
  real *hDiffOrder;  // Hyperdiffusion order.
  real *hDiff;       // Hyperdiffusion coefficient.
  real *kDiffMin;    // Minimum k at which to apply HD.
  // The following are used by initial conditions. It'd be good to put them elsewhere.
  mint icOp;         // IC option.
  real *initAux;     // Auxiliary parameters for ICs.
  real initA;        // Initial amplitude.
  real noiseA;       // Initial noise amplitude.
};

struct mugy_population {
  mint numSpecies;       // Number of species.
  mint mpiProcs;         // MPI decomposition of the species.
  mint globalSpecOff;    // Offset of first species in this process within the global population.
  mint globalMomOff;     // Offset of first moment in this process within global number of moments.
  struct mugy_species *spec;  // Pointer to array of species.
  mint numMomentsTot;    // Total number of moments across all species.
  struct mugy_fourierArray *momk; // Moments in Fourier space. Possibly multiple copies (e.g. for time stepper).
};

// Return a pointer to the momIdx-th moment of the sIdx-th species in mom/momk.
real* getMoment_real(struct mugy_realGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, real *momIn);
fourier* getMoment_fourier(struct mugy_fourierGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, fourier *momkIn);

#endif
