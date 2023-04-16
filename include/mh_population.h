/* mugy: mh_population.h
 *
 * Population object holding the data and information about
 * each species (including its moments).
 *
 */
#pragma once

#include "mh_data.h"
#include "mh_grid.h"
#include "mh_array.h"

struct mugy_species_pars {
  mint numMoments;   // Number of moments.
  real qCharge;      // Charge.
  real muMass;       // sqrt of the reference ion mass to mass of this species.
  real tau;          // Ratio of reference ion temperature to temperature of this species. 
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

struct mugy_pop {
  mint numSpecies;       // Number of species.
  mint globalSpecOff;    // Offset of first species in this process within the global population.
  mint globalMomOff;     // Offset of first moment in this process within global number of moments.
  struct mugy_species_pars *pars;  // Pointer to array of species parameters.
  mint numMomentsTot;    // Total number of moments across all species.
  struct mugy_array **momk;  // Moments in Fourier space. Possibly multiple copies (e.g. for time stepper).
  struct mugy_array *pbFLRop;  // Poisson bracket FLR operators for each species.
  struct mugy_array *poissonFac;  // Factors multiplying each moment in Poissson equation.
};

struct mugy_population {
  struct mugy_pop local;   // Population local to this MPI process.
  struct mugy_pop global;  // Global population.
  mint mpiProcs;         // MPI decomposition of the species.
};

struct mugy_population *mugy_population_alloc();

// Allocate real-space moment vectors, on host and/or device.
//   grid: grid on which to allocate the vector of moments.
//   pop: population struct mugy_containing the number of species and moments.
//   res: resource on which to allocate (host, device or both).
struct mugy_array *mugy_population_alloc_realMoments(const struct mugy_realGrid grid, const struct mugy_pop pop, enum mugy_resource_mem res);

// Allocate Fourier-space moment vectors, on host and/or device.
//   grid: grid on which to allocate the vector of moments.
//   pop: population struct mugy_containing the number of species and moments.
//   res: resource on which to allocate (host, device or both).
struct mugy_array *mugy_population_alloc_fourierMoments(const struct mugy_fourierGrid grid, const struct mugy_pop pop, enum mugy_resource_mem res);

void mugy_population_alloc_moments(struct mugy_population *pop, struct mugy_grid *grid);

// Return a pointer to the momIdx-th moment of the sIdx-th species in mom/momk.
real* mugy_population_getMoment_real(struct mugy_realGrid grid, struct mugy_pop pop, mint sIdx, mint momIdx, real *momIn);
void* mugy_population_getMoment_fourier(struct mugy_fourierGrid grid, struct mugy_pop pop, mint sIdx, mint momIdx, void *momkIn);

void mugy_population_free(struct mugy_population *pop);

