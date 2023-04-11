/* mugy: mh_ffts_dev.h
 * 
 * Compute FFTs on the device.
 *
 */

#ifndef MUGY_FFTS_DEV
#define MUGY_FFTS_DEV

#if USE_GPU

#include "mh_grid.h"
#include "mh_comms.h"
#include "mh_macros.h"


#ifdef __NVCC__

#include <cufftXt.h>
#include "mh_cufft_utils.h"

#if USE_SINGLE_PRECISION
typedef cufftComplex mugy_cufft_fourier;
#else
typedef cufftDoubleComplex mugy_cufft_fourier;
#endif

typedef cufftHandle mugy_cufft_plan;

// Info needed for a single FFT on the device.
struct mugy_fft_dev {
  mugy_cufft_fourier *kbuf;  // Fourier-space buffer (on device).
  real *rbuf;                // Real-space buffer (on device).
  real normFac;              // Normalization.
  bool forwardNorm;          // Normalize in r2c (true) or c2r (false) FFT.
  mugy_cufft_plan plan_r2c, plan_c2r;  // Plans.
};
#endif
// end if __NVCC__

void fft_init_dev(struct mugy_fft_dev *ffts, struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_comms comms);

void fft_xy_c2r_dev(struct mugy_fft_dev *ffts, real *fOut, void *fkIn);
void fft_xy_r2c_dev(struct mugy_fft_dev *ffts, void *fkOut, real *fIn);

void fft_terminate_dev(struct mugy_fft_dev *ffts);

#endif
// end if USE_GPU

#endif
