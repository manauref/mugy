/* mugy: dydt_dev.cu
 *
 * Compute the time rate of change of the moment vector on the device.
 *
 */
extern "C" {
#include "mh_dydt_dev.h"
#include "mh_utilities.h"
#include "mh_array.h"
}
#include "mh_fourier_dev.h"
#include "mh_utilities_dev.h"

// Starting linear index for each thread.
#define LINIDX0 (blockIdx.x*blockDim.x+threadIdx.x)

MUGY_CU_D static inline void* mugy_population_getMomentk_cu(void *momkIn, mint sIdx, mint momIdx, mint xIdx, mint NekxTot, mint numMoments) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species at the xIdx-th grid lication in momk.
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += numMoments;
  return (mugy_cufourier_t *)momkIn + (momOff+momIdx)*NekxTot+xIdx;
}

MUGY_CU_D static inline void* mugy_array_getk_cu(void *arr, mint linIdx) {
  // Get pointer to the element in the array with linear index 'linIdx'.
  return (mugy_cufourier_t *)arr + linIdx;
}

__global__ void mugy_dydt_phikFLR_cu(mugy_cufourier_t *phikFLR, const mugy_cufourier_t *phik,
  const real *pbFLRop, mint nelem, mint numSpecies, mint numMoments, mint NxTot, mint NxyTot) {

  for (mint s=0; s<numSpecies; s++) {
    // Add contributions to the time rate of change of the m-th moment.
    for (unsigned long linIdx = LINIDX0; linIdx < nelem; linIdx += blockDim.x*gridDim.x) {
      const mugy_cufourier_t *phik_p     = (const mugy_cufourier_t *)mugy_array_getk_cu((void*)phik, linIdx);

      unsigned long lin3 = s*3*NxTot + linIdx;
      const mugy_cufourier_t *phikFLR0_p = (const mugy_cufourier_t *)mugy_array_getk_cu((void *)phikFLR, lin3+0*NxyTot);
      const mugy_cufourier_t *phikFLR1_p = (const mugy_cufourier_t *)mugy_array_getk_cu((void *)phikFLR, lin3+1*NxyTot);
      const mugy_cufourier_t *phikFLR2_p = (const mugy_cufourier_t *)mugy_array_getk_cu((void *)phikFLR, lin3+2*NxyTot);

      // Convert the 3D linIdx to a 2D (perp) linIdx to index FLR ops.
//      mint idx[nDim], idxperp[2];
//      mugy_grid_lin2sub_fourier(idx, linIdx, gridL);
//      for (mint d=0; d<2; d++) idxperp[d] = idx[d];
//      unsigned long linIdxperp = mugy_grid_sub2lin_perp_fourier(idxperp, gridL);
//
//      unsigned long lin3perp = s*3*NxyTot + linIdxperp;
//      const real *flrOp0_p = (const real *)mugy_array_get_cu((void*)pbFLRop, lin3perp+0*NxyTot);
//      const real *flrOp1_p = (const real *)mugy_array_get_cu((void*)pbFLRop, lin3perp+1*NxyTot);
//      const real *flrOp2_p = (const real *)mugy_array_get_cu((void*)pbFLRop, lin3perp+2*NxyTot);
//
//      phikFLR0_p[0] = mugy_cuCmul(mugy_make_cuComplex(flrOp0_p[0],0.), phik_p[0]);
//      phikFLR1_p[0] = mugy_cuCmul(mugy_make_cuComplex(flrOp1_p[0],0.), phik_p[0]);
//      phikFLR2_p[0] = mugy_cuCmul(mugy_make_cuComplex(flrOp2_p[0],0.), phik_p[0]);
    }
  }
}
void mugy_dydt_phikFLR_dev(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid) {
  // Compute <J_0>phik, 0.5*hatLap<J_0>phik, 1+0.5*hatLap+hathatLap)<J_0>phik on device.

  struct mugy_population_species *popL = pop->local;
  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_array *phik    = field->phik;
  struct mugy_array *phikFLR = popL->phikFLR;
  struct mugy_array *pbFLRop = popL->pbFLRop;

  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(phik->nelem, nthreads);

  // WARNING: assume all species have the same numMoments.
//  mugy_dydt_phikFLR_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)phikFLR->dev, (const mugy_cufourier_t *)phik->dev,
//    (const real*)pbFLRop->dev, phik->nelem, popL->numSpecies, popL->pars[0].numMoments, gridL->NxTot);
}

__global__ void mugy_dydt_linear_cu(mugy_cufourier_t *momkDot, const mugy_cufourier_t *momk, const mugy_cufourier_t *phik,
  const mugy_cufourier_t *linOpMom, const mugy_cufourier_t *linOpPhi, mint nelem, mint numSpecies, mint numMoments, mint NxTot) {

  mint sOff = 0;
  for (mint s=0; s<numSpecies; s++) {
    for (mint m=0; m<numMoments; m++) {
      // Add contributions to the time rate of change of the m-th moment.
      for (unsigned long linIdx = LINIDX0; linIdx < nelem; linIdx += blockDim.x*gridDim.x) {
	mugy_cufourier_t *momkDot_p = (mugy_cufourier_t *)mugy_population_getMomentk_cu((void*)momkDot, s, m, linIdx, NxTot, numMoments);

	// Linear operator acting on the n-th moment.
        for (mint n=0; n<numMoments; n++) {
	  const mugy_cufourier_t *momk_p     = (const mugy_cufourier_t *)mugy_population_getMomentk_cu((void*)momk, s, n, linIdx, NxTot, numMoments);
          const mugy_cufourier_t *linOpMom_p = (const mugy_cufourier_t *)mugy_array_getk_cu((void*)linOpMom, sOff+(m*numMoments+n)*NxTot+linIdx);
          momkDot_p[0] = mugy_cuCadd(momkDot_p[0], mugy_cuCmul(linOpMom_p[0],momk_p[0]));
        }

	// Linear operator acting on the potential.
        const mugy_cufourier_t *linOpPhi_p = (const mugy_cufourier_t *)mugy_population_getMomentk_cu((void*)linOpPhi, s, m, linIdx, NxTot, numMoments);
        const mugy_cufourier_t *phik_p     = (const mugy_cufourier_t *)mugy_array_getk_cu((void*)phik, linIdx);
        momkDot_p[0] = mugy_cuCadd(momkDot_p[0], mugy_cuCmul(linOpPhi_p[0],phik_p[0]));
      }
    }
    sOff += numMoments*numMoments*NxTot;
  }
}

void mugy_dydt_linear_dev(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid) {
  // Compute the time rate of change dy/dt, where y is the vector of moments, due to linear terms.

  struct mugy_population_species *popL = pop->local;
  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_array *momkDot = popL->momk[outIdx];
  struct mugy_array *momk    = popL->momk[inIdx];
  struct mugy_array *phik    = field->phik;

  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(phik->nelem, nthreads);

  // WARNING: assume all species have the same numMoments.
  mugy_dydt_linear_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)momkDot->dev, (const mugy_cufourier_t *)momk->dev,
    (const mugy_cufourier_t *)phik->dev, (const mugy_cufourier_t *)popL->linOpMom->dev,
    (const mugy_cufourier_t *)popL->linOpPhi->dev, phik->nelem, popL->numSpecies, popL->pars[0].numMoments, gridL->NxTot);
}

void mugy_dydt_dev(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid, struct mugy_fft *fft) {
  // Compute the time rate of change dy/dt, where y is the vector of moments.

  struct mugy_population_species *popL = pop->local;
  struct mugy_array *momkDot = popL->momk[outIdx];
  mugy_array_zero(momkDot, MUGY_DEVICE_MEM);

  // Apply linear terms.
  mugy_dydt_linear_dev(outIdx, inIdx, time, pop, field, grid);
}
