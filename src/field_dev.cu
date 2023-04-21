/* mugy: field_dev.cu
 *
 * Manipulate the field object on the device.
 *
 */

extern "C" {
#include "mh_field_dev.h"
#include "mh_utilities.h"
}
#include "mh_fourier_dev.h"
#include "mh_utilities_dev.h"

// Starting linear index for each thread.
#define LINIDX0 (blockIdx.x*blockDim.x+threadIdx.x)

MUGY_CU_D void* getMomentk(void *momkIn, mint sIdx, mint momIdx, mint xIdx, mint NxTot, mint numMoments) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species at the xIdx-th grid lication in momk.
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += numMoments;
  return (mugy_cufourier_t *)momkIn + (momOff+momIdx)*NxTot+xIdx;
}

__global__ void mugy_field_poisson_solve_cu(mugy_cufourier_t *phik, const mugy_cufourier_t *momk,
  const mugy_cufourier_t *poissonFac, mint nelem, mint numSpecies, mint numMoments, mint NxTot)
{
  for (unsigned long linIdx = LINIDX0; linIdx < nelem; linIdx += blockDim.x*gridDim.x) {
    mugy_cufourier_t *phik_p = &phik[linIdx];
    phik_p[0] = mugy_make_cuComplex(0., 0.);
    for (mint s=0; s<numMoments; s++) {
      for (mint m=0; m<numMoments; m++) {
        const mugy_cufourier_t *poissonFac_p = (const mugy_cufourier_t *)getMomentk((void*)poissonFac, s, m, linIdx, NxTot, numMoments);
        const mugy_cufourier_t *momk_p       = (const mugy_cufourier_t *)getMomentk((void*)momk      , s, m, linIdx, NxTot, numMoments);

        phik_p[0] = mugy_cuCadd(phik_p[0], mugy_cuCmul(poissonFac_p[0],momk_p[0]));
      }
    }
  }
}

void mugy_field_poisson_solve_dev(struct mugy_field *field, struct mugy_population *pop, struct mugy_grid *grid, mint tstepIdx) {
  // Solve the Poisson equation to obtain phik using the charge density
  // in the time-stepping index 'tstepIdx'

  struct mugy_population_species *popL = pop->local;
  struct mugy_grid_basic *gridL = grid->local->fourier;
  struct mugy_array      *phik  = field->phik;
  struct mugy_array      *momk  = popL->momk[tstepIdx];
  
  mint nthreads = DEFAULT_NUM_THREADS_DEV;
  mint nblocks  = mugy_div_up_mint(phik->nelem, nthreads);

  // WARNING: assume all species have the same numMoments.
  mugy_field_poisson_solve_cu<<<nblocks, nthreads>>>((mugy_cufourier_t *)phik->dev,
    (const mugy_cufourier_t *)momk->dev, (const mugy_cufourier_t *)popL->poissonFac->dev,
    phik->nelem, popL->numSpecies, popL->pars[0].numMoments, gridL->NxTot);
}
