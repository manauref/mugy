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
#include <stdbool.h>
#include <stdlib.h>  // for malloc.

// Info needed for a single FFT on the host.
struct mugy_fft_ho {
  mugy_fftw_fourier *kbuf;  // Fourier-space buffer.
  real *rbuf;               // Real-space buffer.
  real normFac;             // Normalization.
  bool forwardNorm;         // Normalize in r2c (true) or c2r (false) FFT.
  mugy_fftw_plan plan_r2c, plan_c2r;  // Plans.
};

struct mugy_fft_fam_ho {
  struct mugy_fft_ho *xy;    // 2D FFT of x-y planes.
  struct mugy_fft_ho *xy_a;  // 2D FFT of (aliased) x-y planes.
};

// Info needed for all FFTs in mugy.
struct mugy_ffts {
  struct mugy_fft_fam_ho *ho; 
  struct mugy_fft_fam_dev *dev;
};

struct mugy_ffts *fft_init(struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_comms comms) {

  // Allocate space for the FFT manager.
  struct mugy_ffts *ffts = (struct mugy_ffts *) malloc(sizeof(struct mugy_ffts));

#if USE_GPU
  // Initialize device FFTs.
  ffts->dev = mugy_fft_init_dev(gridG, gridL, comms);
#endif

  // Initialize host FFTs.
  mugy_fftw_mpi_init();

  ffts->ho = (struct mugy_fft_fam_ho *) malloc(sizeof(struct mugy_fft_fam_ho));

  // ....... Setup for 2D FFTs of a single field ......... //
  ffts->ho->xy = (struct mugy_fft_ho *) malloc(sizeof(struct mugy_fft_ho));
  struct mugy_fft_ho *cfft = ffts->ho->xy;  // Temp pointer for convenience.
  struct mugy_comms_sub *scomm = &comms.sub2d[0];  // Temp pointer for convenience.
  // Get local data size.
  const mint fftDim = 2;
  ptrdiff_t fftSizek[fftDim], fftNum, blockSizek0;
  fftSizek[0] = gridG.fG.Nekx[0];
  fftSizek[1] = gridG.fG.Nekx[1];
  fftNum      = gridL.fG.Nekx[2];
  blockSizek0 = gridL.fG.Nekx[0];
  ptrdiff_t alloc_local, local_Nekx0, local_kx0_start;
  alloc_local = mugy_fftw_mpi_local_size_many(fftDim, fftSizek, fftNum, blockSizek0, scomm->comm, &local_Nekx0, &local_kx0_start);

  // Allocate buffers.
  cfft->kbuf = mugy_fftw_alloc_complex(alloc_local);
  cfft->rbuf = mugy_fftw_alloc_real(2*alloc_local);

  // Create plans.
  ptrdiff_t fftSize[fftDim], blockSize0in, blockSize0out;
  fftSize[0] = gridG.fG.dual.Nx[0];
  fftSize[1] = gridG.fG.dual.Nx[1];
  blockSize0in  = gridL.fG.dual.Nx[0];
  blockSize0out = blockSizek0;
  cfft->plan_r2c = mugy_fftw_mpi_plan_many_dft_r2c(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                      cfft->rbuf, cfft->kbuf, scomm->comm, FFTW_ESTIMATE);
  blockSize0in  = blockSizek0;
  blockSize0out = gridL.fG.dual.Nx[0];
  cfft->plan_c2r = mugy_fftw_mpi_plan_many_dft_c2r(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                      cfft->kbuf, cfft->rbuf, scomm->comm, FFTW_ESTIMATE);

  cfft->normFac = 1./((real)gridG.fG.dual.NxyTot);
  cfft->forwardNorm = false;  // This FFT is only used for ICs given in real-space.
  // ....... End setup for 2D FFTs of a single field ......... //

  return ffts;
}

void fft_xy_c2r(struct mugy_ffts *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum resource_comp res) {
#if USE_GPU
//  if (res == deviceComp)
//    return mugy_fft_xy_c2r_dev(ffts->dev, fOut, fkIn);
#endif
  
  struct mugy_fft_ho *cfft = ffts->ho->xy;  // Temp pointer for convenience.

  // Copy data into buffer.
  memcpy_fourier(cfft->kbuf, fkIn->ho, fkIn->nelem, host2host);

  // Inverse FFT.
  mugy_fftw_mpi_execute_dft_c2r(cfft->plan_c2r, cfft->kbuf, cfft->rbuf);

  // Copy data from buffer.
  memcpy_real(fOut->ho, cfft->rbuf, fOut->nelem, host2host);

  // Apply the nonunitary normalization.
  if (!cfft->forwardNorm)
    mugy_array_scale(fOut, cfft->normFac, hostComp);
}

void fft_xy_r2c(struct mugy_ffts *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum resource_comp res) {
#if USE_GPU
//  if (res == deviceComp)
//    return mugy_fft_xy_r2c_dev(ffts->dev, fkOut, fIn);
#endif
  
  struct mugy_fft_ho *cfft = ffts->ho->xy;  // Temp pointer for convenience.

  // Copy data into buffer.
  memcpy_real(cfft->rbuf, fIn->ho, fIn->nelem, host2host);

  // Forward FFT.
  mugy_fftw_mpi_execute_dft_r2c(cfft->plan_r2c, cfft->rbuf, cfft->kbuf);

  // Copy data from buffer.
  memcpy_fourier(fkOut->ho, cfft->kbuf, fkOut->nelem, host2host);

  // Apply the nonunitary normalization.
  if (cfft->forwardNorm)
    mugy_array_scale(fkOut, cfft->normFac, hostComp);
}

void fft_terminate(struct mugy_ffts *ffts) {
#if USE_GPU
  mugy_fft_terminate_dev(ffts->dev);  // Free memory for device FFTs.
#endif

  // Deallocate objects for xy FFT of a single field.
  struct mugy_fft_ho *cfft = ffts->ho->xy;  // Temp pointer for convenience.
  mugy_fftw_free(cfft->kbuf);
  mugy_fftw_free(cfft->rbuf);
  mugy_fftw_destroy_plan(cfft->plan_c2r);
  mugy_fftw_destroy_plan(cfft->plan_r2c);
  free(cfft);

  free(ffts->ho);
  free(ffts);
}
