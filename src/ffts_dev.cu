/* mugy: ffts_dev.cu
 *
 * Device FFT methods.
 *
 */

extern "C" {
#include "mh_ffts_dev.h"
}
#include "mh_ffts_dev_priv.h"

struct mugy_fft_fam_dev *mugy_fft_init_dev(struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_comms comms) {

  struct mugy_fft_fam_dev *ffts = (struct mugy_fft_fam_dev *) malloc(sizeof(struct mugy_fft_fam_dev));

  // ....... Setup for 2D FFTs of a single field ......... //
  ffts->xy = (struct mugy_fft_dev *) malloc(sizeof(struct mugy_fft_dev));
  struct mugy_fft_dev *cfft = ffts->xy;  // Temp pointer for convenience.

  CUDA_RT_CALL(cudaStreamCreateWithFlags(&cfft->stream, cudaStreamNonBlocking));

  size_t workSizes_r2c[1], workSizes_c2r[1];
  const mint fftDim = 2;
  long long int fftSize[fftDim], fftNum;
  long long int inembed[nDim], istride, idist;
  long long int onembed[nDim], ostride, odist;
  fftSize[0] = (long long int) gridG.fG.dual.Nx[0];
  fftSize[1] = (long long int) gridG.fG.dual.Nx[1];
  fftNum     = gridL.fG.Nekx[2];
  istride    = 1;
  ostride    = 1;

  inembed[0] = (long long int) gridL.fG.Nekx[2];
  inembed[1] = (long long int) gridL.fG.dual.Nx[0];
  inembed[2] = (long long int) gridL.fG.dual.Nx[1];
  idist      = (long long int) gridL.fG.dual.NxyTot;
  onembed[0] = (long long int) gridL.fG.Nekx[2];
  onembed[1] = (long long int) gridL.fG.Nekx[0];
  onembed[2] = (long long int) gridL.fG.Nekx[1];
  odist      = (long long int) gridL.fG.NekxyTot;
  CUFFT_CALL(cufftCreate(&cfft->plan_r2c));
  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_r2c, fftDim, fftSize,
    inembed, istride, idist, mugy_cufft_real_enum,
    onembed, ostride, odist, mugy_cufft_fourier_enum,
    fftNum, workSizes_r2c, mugy_cufft_executiontype));
  CUFFT_CALL(cufftSetStream(cfft->plan_r2c, cfft->stream));

  inembed[0] = (long long int) gridL.fG.Nekx[2];
  inembed[1] = (long long int) gridL.fG.Nekx[0];
  inembed[2] = (long long int) gridL.fG.Nekx[1];
  idist      = (long long int) gridL.fG.NekxyTot;
  onembed[0] = (long long int) gridL.fG.Nekx[2];
  onembed[1] = (long long int) gridL.fG.dual.Nx[0];
  onembed[2] = (long long int) gridL.fG.dual.Nx[1];
  odist      = (long long int) gridL.fG.dual.NxyTot;
  CUFFT_CALL(cufftCreate(&cfft->plan_c2r));
  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_c2r, fftDim, fftSize,
    inembed, istride, idist, mugy_cufft_fourier_enum,
    onembed, ostride, odist, mugy_cufft_real_enum,
    fftNum, workSizes_c2r, mugy_cufft_executiontype));
  CUFFT_CALL(cufftSetStream(cfft->plan_c2r, cfft->stream));

  // Allocate buffers.
  // Switch these mallocs to mugy functions.
  CUDA_RT_CALL(cudaMalloc(&cfft->rbuf, sizeof(real) * gridL.fG.dual.NxTot));
  CUDA_RT_CALL(cudaMalloc(&cfft->kbuf, sizeof(mugy_cufft_fourier) * gridL.fG.NekxTot));

  cfft->normFac = 1./((real)gridG.fG.dual.NxyTot);
  cfft->forwardNorm = false;  // This FFT is only used for ICs given in real-space.
  // ....... End setup for 2D FFTs of a single field ......... //

  return ffts;
}

void mugy_fft_xy_c2r_dev(struct mugy_fft_fam_dev *ffts, struct mugy_array *fOut, struct mugy_array *fkIn) {

  struct mugy_fft_dev *cfft = ffts->xy;  // Temp pointer for convenience.

  // Copy data into buffer.
  mugy_memcpy(cfft->kbuf, fkIn->dev, fkIn->nelemsz, device2device);

  // Inverse FFT.
  cufftXtExec(cfft->plan_c2r, cfft->kbuf, cfft->rbuf, 0);

  // Copy data from buffer.
  mugy_memcpy(fOut->dev, cfft->rbuf, fOut->nelemsz, device2device);

  // Apply the nonunitary normalization.
  if (!cfft->forwardNorm)
    mugy_array_scale(fOut, cfft->normFac, deviceComp);
  
}

void mugy_fft_xy_r2c_dev(struct mugy_fft_fam_dev *ffts, struct mugy_array *fkOut, struct mugy_array *fIn) {

  struct mugy_fft_dev *cfft = ffts->xy;  // Temp pointer for convenience.

  // Copy data into buffer.
  mugy_memcpy(cfft->rbuf, fIn->dev, fIn->nelemsz, device2device);

  // Forward FFT.
  cufftXtExec(cfft->plan_r2c, cfft->rbuf, cfft->kbuf, 0);

  // Copy data from buffer.
  mugy_memcpy(fkOut->dev, cfft->kbuf, fkOut->nelemsz, device2device);

  // Apply the nonunitary normalization.
  if (cfft->forwardNorm)
    mugy_array_scale(fkOut, cfft->normFac, deviceComp);
  
}

void mugy_fft_terminate_dev(struct mugy_fft_fam_dev *ffts) {

  struct mugy_fft_dev *cfft = ffts->xy;  // Temp pointer for convenience.
  CUFFT_CALL(cufftDestroy(cfft->plan_r2c));
  CUFFT_CALL(cufftDestroy(cfft->plan_c2r));
  CUDA_RT_CALL(cudaStreamDestroy(cfft->stream));
  CUDA_RT_CALL(cudaFree(cfft->rbuf));
  CUDA_RT_CALL(cudaFree(cfft->kbuf));
  free(cfft);

  free(ffts);
}
