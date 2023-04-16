/* mugy: ffts_dev.cu
 *
 * Device FFT methods.
 *
 */

extern "C" {
#include "mh_ffts.h"
#include "mh_ffts_dev.h"
#include "mh_utilities.h"
}
#include "mh_ffts_dev_priv.h"

struct mugy_fft_fam_dev *mugy_fft_init_dev(struct mugy_grid *grid, struct mugy_pop popL, struct mugy_comms *comms) {

  struct mugy_fft_fam_dev *ffts = (struct mugy_fft_fam_dev *) malloc(sizeof(struct mugy_fft_fam_dev));

  struct mugy_fft_dev *cfft;  // Temp pointer for convenience.

  // ....... Setup for 2D FFTs of a single field ......... //
  ffts->xy = (struct mugy_fft_dev *) malloc(sizeof(struct mugy_fft_dev));
  cfft = ffts->xy;

  CUDA_RT_CALL(cudaStreamCreateWithFlags(&cfft->stream, cudaStreamNonBlocking));

  size_t workSizes_r2c[1], workSizes_c2r[1];
  const mint fftDim = 2;
  long long int fftSize[fftDim], fftNum;
  long long int inembed[fftDim], istride, idist;
  long long int onembed[fftDim], ostride, odist;
  fftSize[0] = grid->global.deal.dual.Nx[0];
  fftSize[1] = grid->global.deal.dual.Nx[1];
  fftNum     = grid->local.deal.Nekx[2];
  istride    = 1;
  ostride    = 1;

  inembed[0] = grid->local.deal.dual.Nx[0];
  inembed[1] = grid->local.deal.dual.Nx[1];
  idist      = grid->local.deal.dual.NxyTot;
  onembed[0] = grid->local.deal.Nekx[0];
  onembed[1] = grid->local.deal.Nekx[1];
  odist      = grid->local.deal.NekxyTot;
  CUFFT_CALL(cufftCreate(&cfft->plan_r2c));
  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_r2c, fftDim, fftSize,
    inembed, istride, idist, MUGY_CUFFT_REAL,
    onembed, ostride, odist, MUGY_CUFFT_FOURIER,
    fftNum, workSizes_r2c, MUGY_CUFFT_EXEC_FOURIER));
  CUFFT_CALL(cufftSetStream(cfft->plan_r2c, cfft->stream));

  inembed[0] = grid->local.deal.Nekx[0];
  inembed[1] = grid->local.deal.Nekx[1];
  idist      = grid->local.deal.NekxyTot;
  onembed[0] = grid->local.deal.dual.Nx[0];
  onembed[1] = grid->local.deal.dual.Nx[1];
  odist      = grid->local.deal.dual.NxyTot;
  CUFFT_CALL(cufftCreate(&cfft->plan_c2r));
  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_c2r, fftDim, fftSize,
    inembed, istride, idist, MUGY_CUFFT_FOURIER,
    onembed, ostride, odist, MUGY_CUFFT_REAL,
    fftNum, workSizes_c2r, MUGY_CUFFT_EXEC_FOURIER));
  CUFFT_CALL(cufftSetStream(cfft->plan_c2r, cfft->stream));

  // Allocate buffers.
  // Switch these mallocs to mugy functions.
  CUDA_RT_CALL(cudaMalloc(&cfft->rbuf, grid->local.deal.dual.NxTot*sizeof(real)));
  CUDA_RT_CALL(cudaMalloc(&cfft->kbuf, grid->local.deal.NekxTot*sizeof(mugy_cufft_fourier_t)));

  cfft->normFac = 1./((real)grid->global.deal.dual.NxyTot);
  cfft->forwardNorm = false;  // This FFT is only used for ICs given in real-space.

  // ....... End setup for 2D FFTs of a single field ......... //

  // ....... Setup for xy FFTs of all moments ......... //
  ffts->mom_xy = (struct mugy_fft_dev *) malloc(sizeof(struct mugy_fft_dev));
  cfft = ffts->mom_xy;

  CUDA_RT_CALL(cudaStreamCreateWithFlags(&cfft->stream, cudaStreamNonBlocking));

  fftNum     = popL.numMomentsTot*grid->local.deal.Nekx[2];

  inembed[0] = grid->local.deal.dual.Nx[0];
  inembed[1] = grid->local.deal.dual.Nx[1];
  idist      = grid->local.deal.dual.NxyTot;
  onembed[0] = grid->local.deal.Nekx[0];
  onembed[1] = grid->local.deal.Nekx[1];
  odist      = grid->local.deal.NekxyTot;
  CUFFT_CALL(cufftCreate(&cfft->plan_r2c));
  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_r2c, fftDim, fftSize,
    inembed, istride, idist, MUGY_CUFFT_REAL,
    onembed, ostride, odist, MUGY_CUFFT_FOURIER,
    fftNum, workSizes_r2c, MUGY_CUFFT_EXEC_FOURIER));
  CUFFT_CALL(cufftSetStream(cfft->plan_r2c, cfft->stream));

  inembed[0] = grid->local.deal.Nekx[0];
  inembed[1] = grid->local.deal.Nekx[1];
  idist      = grid->local.deal.NekxyTot;
  onembed[0] = grid->local.deal.dual.Nx[0];
  onembed[1] = grid->local.deal.dual.Nx[1];
  odist      = grid->local.deal.dual.NxyTot;
  CUFFT_CALL(cufftCreate(&cfft->plan_c2r));
  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_c2r, fftDim, fftSize,
    inembed, istride, idist, MUGY_CUFFT_FOURIER,
    onembed, ostride, odist, MUGY_CUFFT_REAL,
    fftNum, workSizes_c2r, MUGY_CUFFT_EXEC_FOURIER));
  CUFFT_CALL(cufftSetStream(cfft->plan_c2r, cfft->stream));

  // Allocate buffers.
  // Switch these mallocs to mugy functions.
  CUDA_RT_CALL(cudaMalloc(&cfft->rbuf, popL.numMomentsTot*grid->local.deal.dual.NxTot*sizeof(real)));
  CUDA_RT_CALL(cudaMalloc(&cfft->kbuf, popL.numMomentsTot*grid->local.deal.NekxTot*sizeof(mugy_cufft_fourier_t)));

  cfft->normFac = 1./((real)grid->global.deal.dual.NxyTot);
  cfft->forwardNorm = false;  // This FFT is only used for ICs given in real-space.

  // ....... End setup for xy FFTs of all moments ......... //

  return ffts;
}

void mugy_fft_c2r_dev(struct mugy_fft_fam_dev *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum mugy_fft_type ttype) {

  struct mugy_fft_dev *cfft;  // Temp pointer for convenience.
  if (ttype == mugy_fft_xy)
    cfft = ffts->xy;
  else if (ttype == mugy_fft_mom_xy)
    cfft = ffts->mom_xy;
  else
    abortSimulation(" mugy_fft_c2r_dev: fft type not supported! Terminating...\n");

  // Copy data into buffer.
  mugy_memcpy(cfft->kbuf, fkIn->dev, fkIn->nelemsz, MUGY_DEVICE2DEVICE);

  // Inverse FFT.
  CUFFT_CALL(cufftXtExec(cfft->plan_c2r, cfft->kbuf, cfft->rbuf, CUFFT_INVERSE));

  CUDA_RT_CALL(cudaStreamSynchronize(cfft->stream));

  // Copy data from buffer.
  mugy_memcpy(fOut->dev, cfft->rbuf, fOut->nelemsz, MUGY_DEVICE2DEVICE);

  // Apply the nonunitary normalization.
  if (!cfft->forwardNorm)
    mugy_array_scale(fOut, cfft->normFac, MUGY_DEVICE_CALC);
  
}

void mugy_fft_r2c_dev(struct mugy_fft_fam_dev *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum mugy_fft_type ttype) {

  struct mugy_fft_dev *cfft;  // Temp pointer for convenience.
  if (ttype == mugy_fft_xy)
    cfft = ffts->xy;
  else if (ttype == mugy_fft_mom_xy)
    cfft = ffts->mom_xy;
  else
    abortSimulation(" mugy_fft_r2c_dev: fft type not supported! Terminating...\n");

  // Copy data into buffer.
  mugy_memcpy(cfft->rbuf, fIn->dev, fIn->nelemsz, MUGY_DEVICE2DEVICE);

  // Forward FFT.
  CUFFT_CALL(cufftXtExec(cfft->plan_r2c, cfft->rbuf, cfft->kbuf, CUFFT_FORWARD));

  CUDA_RT_CALL(cudaStreamSynchronize(cfft->stream));

  // Copy data from buffer.
  mugy_memcpy(fkOut->dev, cfft->kbuf, fkOut->nelemsz, MUGY_DEVICE2DEVICE);

  // Apply the nonunitary normalization.
  if (cfft->forwardNorm)
    mugy_array_scale(fkOut, cfft->normFac, MUGY_DEVICE_CALC);
  
}

void mugy_fft_terminate_dev(struct mugy_fft_fam_dev *ffts) {

  struct mugy_fft_dev *cfft = ffts->xy;  // Temp pointer for convenience.

  cfft = ffts->xy;
  CUDA_RT_CALL(cudaStreamDestroy(cfft->stream));
  CUFFT_CALL(cufftDestroy(cfft->plan_r2c));
  CUFFT_CALL(cufftDestroy(cfft->plan_c2r));
  CUDA_RT_CALL(cudaFree(cfft->rbuf));
  CUDA_RT_CALL(cudaFree(cfft->kbuf));
  free(cfft);

  cfft = ffts->mom_xy;
  CUDA_RT_CALL(cudaStreamDestroy(cfft->stream));
  CUFFT_CALL(cufftDestroy(cfft->plan_r2c));
  CUFFT_CALL(cufftDestroy(cfft->plan_c2r));
  CUDA_RT_CALL(cudaFree(cfft->rbuf));
  CUDA_RT_CALL(cudaFree(cfft->kbuf));
  free(cfft);

  free(ffts);
}
