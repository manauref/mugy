/* mugy: ffts.c
 *
 * Module with FFT operators. We do host FFTs with FFTW
 * and device FFTs with cuFFT. We also anticipate 3 modes
 * depending on MPI decomposition (recall data is organized
 * as (z,x,y) in row-major order):
 *   FFTW:
 *     - serial or MPI decomposed along x.
 *     - MPI decomposed along y.
 *     - MPI decomposed alog x and y.
 *   cuFFT (don't yet know the data layout):
 *     - single GPU 
 *     - single node
 *     - multi node
 *
 */

#include "mh_ffts.h"
#include "mh_fftw_wrap.h"
#include "mh_ffts_dev.h"
#include "mh_data.h"
#include "mh_utilities.h"
#include <stdbool.h>
#include <stdlib.h>  // for malloc.

// Info needed for a single FFT on the host.
struct mugy_fft_ho {
  mugy_fftw_fourier_t *kbuf;  // Fourier-space buffer.
  real *rbuf;                 // Real-space buffer.
  real normFac;               // Normalization.
  bool forwardNorm;           // Normalize in r2c (true) or c2r (false) FFT.
  mugy_fftw_plan_t plan_r2c, plan_c2r;  // Plans.
};

struct mugy_fft_fam_ho {
  // FFTs on dealised grids.
  struct mugy_fft_ho *xy;      // xy FFT of single 3D array.
  struct mugy_fft_ho *mom_xy;  // xy FFT of all moments.
  // FFTs on aliased grids.
  struct mugy_fft_ho *xy_a;      // xy FFT of single 3D array.
  struct mugy_fft_ho *mom_xy_a;  // xy FFT of all moments.
};

// Info needed for all FFTs in mugy.
struct mugy_ffts {
  struct mugy_fft_fam_ho *ho; 
  struct mugy_fft_fam_dev *dev;
};

struct mugy_ffts *mugy_fft_init(struct mugy_grid *grid, struct mugy_pop popL, struct mugy_comms *comms) {

  // Allocate space for the FFT manager.
  struct mugy_ffts *ffts = (struct mugy_ffts *) malloc(sizeof(struct mugy_ffts));

#if USE_GPU
  // Initialize device FFTs.
  ffts->dev = mugy_fft_init_dev(grid, popL, comms);
#endif

  // Initialize host FFTs.
  mugy_fftw_mpi_init();

  ffts->ho = (struct mugy_fft_fam_ho *) malloc(sizeof(struct mugy_fft_fam_ho));

  struct mugy_fft_ho *cfft;  // Temp pointer for convenience.
  struct mugy_comms_sub *scomm;  // Temp pointer for convenience.

  // ....... Setup for xy FFTs of a single 3D array ......... //
  ffts->ho->xy = (struct mugy_fft_ho *) malloc(sizeof(struct mugy_fft_ho));
  cfft  = ffts->ho->xy;
  scomm = &comms->sub2d[0];
  // Get local data size.
  const mint fftDim = 2;
  ptrdiff_t fftSizek[fftDim], fftNum, blockSizek0;
  fftSizek[0] = grid->global->fourier->Nx[0];
  fftSizek[1] = grid->global->fourier->Nx[1];
  fftNum      = grid->local->fourier->Nx[2];
  blockSizek0 = grid->local->fourier->Nx[0];
  ptrdiff_t alloc_local, local_Nx0, local_kx0_start;
  alloc_local = mugy_fftw_mpi_local_size_many(fftDim, fftSizek, fftNum, blockSizek0, scomm->comm, &local_Nx0, &local_kx0_start);

  // Allocate buffers.
  cfft->kbuf = mugy_fftw_alloc_complex(alloc_local);
  cfft->rbuf = mugy_fftw_alloc_real(2*alloc_local);

  // Create plans.
  ptrdiff_t fftSize[fftDim], blockSize0in, blockSize0out;
  fftSize[0] = grid->global->real->Nx[0];
  fftSize[1] = grid->global->real->Nx[1];
  blockSize0in  = grid->local->real->Nx[0];
  blockSize0out = blockSizek0;
  cfft->plan_r2c = mugy_fftw_mpi_plan_many_dft_r2c(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                   cfft->rbuf, cfft->kbuf, scomm->comm, FFTW_ESTIMATE);
  blockSize0in  = blockSizek0;
  blockSize0out = grid->local->real->Nx[0];
  cfft->plan_c2r = mugy_fftw_mpi_plan_many_dft_c2r(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                   cfft->kbuf, cfft->rbuf, scomm->comm, FFTW_ESTIMATE);

  cfft->normFac = 1./((real)grid->global->real->NxyTot);
  cfft->forwardNorm = false;  // This FFT is only used for ICs given in real-space.
  // ....... End setup for xy FFTs of a single 3D array ......... //

  // ....... Setup for xy FFTs of all moments ......... //
  ffts->ho->mom_xy = (struct mugy_fft_ho *) malloc(sizeof(struct mugy_fft_ho));
  cfft = ffts->ho->mom_xy;
  fftNum      = popL.numMomentsTot * grid->local->fourier->Nx[2];
  alloc_local = mugy_fftw_mpi_local_size_many(fftDim, fftSizek, fftNum, blockSizek0, scomm->comm, &local_Nx0, &local_kx0_start);

  // Allocate buffers.
  cfft->kbuf = mugy_fftw_alloc_complex(alloc_local);
  cfft->rbuf = mugy_fftw_alloc_real(2*alloc_local);

  // Create plans.
  blockSize0in  = grid->local->real->Nx[0];
  blockSize0out = blockSizek0;
  cfft->plan_r2c = mugy_fftw_mpi_plan_many_dft_r2c(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                   cfft->rbuf, cfft->kbuf, scomm->comm, FFTW_ESTIMATE);
  blockSize0in  = blockSizek0;
  blockSize0out = grid->local->real->Nx[0];
  cfft->plan_c2r = mugy_fftw_mpi_plan_many_dft_c2r(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                   cfft->kbuf, cfft->rbuf, scomm->comm, FFTW_ESTIMATE);

  cfft->normFac = 1./((real)grid->global->real->NxyTot);
  cfft->forwardNorm = false;  // This FFT is only used for ICs given in real-space.
  // ....... End setup for xy FFTs of all moments ......... //

  return ffts;
}

void mugy_fft_c2r(struct mugy_ffts *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum mugy_fft_type ttype, enum mugy_resource_calc res) {
#if USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_fft_c2r_dev(ffts->dev, fOut, fkIn, ttype);
#endif
  
  struct mugy_fft_ho *cfft;  // Temp pointer for convenience.
  if (ttype == mugy_fft_xy)
    cfft = ffts->ho->xy;
  else if (ttype == mugy_fft_mom_xy)
    cfft = ffts->ho->mom_xy;
  else
    abortSimulation(" mugy_fft_c2r: fft type not supported! Terminating...\n");

  // Copy data into buffer.
  mugy_memcpy(cfft->kbuf, fkIn->ho, fkIn->nelemsz, MUGY_HOST2HOST);

  // Inverse FFT.
  mugy_fftw_mpi_execute_dft_c2r(cfft->plan_c2r, cfft->kbuf, cfft->rbuf);

  // Copy data from buffer.
  mugy_memcpy(fOut->ho, cfft->rbuf, fOut->nelemsz, MUGY_HOST2HOST);

  // Apply the nonunitary normalization.
  if (!cfft->forwardNorm)
    mugy_array_scale(fOut, cfft->normFac, MUGY_HOST_CALC);
}

void mugy_fft_r2c(struct mugy_ffts *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum mugy_fft_type ttype, enum mugy_resource_calc res) {
#if USE_GPU
  if (res == MUGY_DEVICE_CALC)
    return mugy_fft_r2c_dev(ffts->dev, fkOut, fIn, ttype);
#endif
  
  struct mugy_fft_ho *cfft;  // Temp pointer for convenience.
  if (ttype == mugy_fft_xy)
    cfft = ffts->ho->xy;
  else if (ttype == mugy_fft_mom_xy)
    cfft = ffts->ho->mom_xy;
  else
    abortSimulation(" mugy_fft_r2c: fft type not supported! Terminating...\n");

  // Copy data into buffer.
  mugy_memcpy(cfft->rbuf, fIn->ho, fIn->nelemsz, MUGY_HOST2HOST);

  // Forward FFT.
  mugy_fftw_mpi_execute_dft_r2c(cfft->plan_r2c, cfft->rbuf, cfft->kbuf);

  // Copy data from buffer.
  mugy_memcpy(fkOut->ho, cfft->kbuf, fkOut->nelemsz, MUGY_HOST2HOST);

  // Apply the nonunitary normalization.
  if (cfft->forwardNorm)
    mugy_array_scale(fkOut, cfft->normFac, MUGY_HOST_CALC);
}

void mugy_fft_terminate(struct mugy_ffts *ffts) {
#if USE_GPU
  mugy_fft_terminate_dev(ffts->dev);  // Free memory for device FFTs.
#endif

  struct mugy_fft_ho *cfft;  // Temp pointer for convenience.

  // Deallocate objects for xy FFT of a single field.
  cfft = ffts->ho->xy;
  mugy_fftw_free(cfft->kbuf);
  mugy_fftw_free(cfft->rbuf);
  mugy_fftw_destroy_plan(cfft->plan_c2r);
  mugy_fftw_destroy_plan(cfft->plan_r2c);
  free(cfft);

  // Deallocate objects for xy FFT of all moments.
  cfft = ffts->ho->mom_xy;
  mugy_fftw_free(cfft->kbuf);
  mugy_fftw_free(cfft->rbuf);
  mugy_fftw_destroy_plan(cfft->plan_c2r);
  mugy_fftw_destroy_plan(cfft->plan_r2c);
  free(cfft);

  free(ffts->ho);
  free(ffts);
}
