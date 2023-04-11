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

//  const mint fftDim = 2;
//  long long int fftSize[fftDim], fftNum;
//  fftSize[0] = (long long int) gridG.fG.dual.Nx[0];
//  fftSize[1] = (long long int) gridG.fG.dual.Nx[1];
//  fftNum     = gridL.fG.Nekx[2];
//  long long int inembed[nDim], istride, idist;
//  inembed[0] = (long long int) gridL.fG.Nekx[2];
//  inembed[1] = (long long int) gridL.fG.dual.Nx[0];
//  inembed[2] = (long long int) gridL.fG.dual.Nx[1];
//  istride    = 1;
//  idist      = (long long int) gridL.fG.dual.NxyTot;
//  long long int onembed[nDim], ostride, odist;
//  onembed[0] = (long long int) gridL.fG.Nekx[2];
//  onembed[1] = (long long int) gridL.fG.Nekx[0];
//  onembed[2] = (long long int) gridL.fG.Nekx[1];
//  ostride    = 1;
//  odist      = (long long int) gridL.fG.NekxyTot;
//
//  size_t workSizes[1];
//
//  CUFFT_CALL(cufftXtMakePlanMany(cfft->plan_r2c, fftDim, fftSize,
//    inembed, istride, idist, mugy_cufft_real_enum,
//    onembed, ostride, odist, mugy_cufft_fourier_enum,
//    fftNum, workSizes, mugy_cufft_executiontype));
  CUFFT_CALL(cufftPlan2d(&cfft->plan_r2c, gridL.fG.dual.Nx[0], gridL.fG.dual.Nx[1], CUFFT_D2Z));

  return ffts;
}

void mugy_fft_xy_c2r_dev(struct mugy_fft_fam_dev *ffts, real *fOut, void *fkIn) {
  
}

void mugy_fft_xy_r2c_dev(struct mugy_fft_fam_dev *ffts, void *fkOut, real *fIn) {
}

void mugy_fft_terminate_dev(struct mugy_fft_fam_dev *ffts) {

  struct mugy_fft_dev *cfft = ffts->xy;  // Temp pointer for convenience.
  CUFFT_CALL(cufftDestroy(cfft->plan_r2c));
  free(cfft);

  free(ffts);
}
