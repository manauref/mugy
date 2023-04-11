/* mugy: mh_ffts_dev_priv.h
 * 
 * Private declarations for device FFTs.
 *
 */
#pragma once

#include "mh_macros.h"
#include <cufftXt.h>
#include "mh_cufft_utils.h"

#if USE_SINGLE_PRECISION
typedef cufftComplex mugy_cufft_fourier;
const cudaDataType mugy_cufft_real_enum = CUDA_R_32F;
const cudaDataType mugy_cufft_fourier_enum = CUDA_C_32F;
const cudaDataType mugy_cufft_executiontype = CUDA_C_32F;
#else
typedef cufftDoubleComplex mugy_cufft_fourier;
const cudaDataType mugy_cufft_real_enum = CUDA_R_64F;
const cudaDataType mugy_cufft_fourier_enum = CUDA_C_64F;
const cudaDataType mugy_cufft_executiontype = CUDA_C_64F;
#endif

typedef cufftHandle mugy_cufft_plan;
typedef cudaStream_t mugy_custream;

// Info needed for a single FFT on the device.
struct mugy_fft_dev {
  mugy_cufft_fourier *kbuf;  // Fourier-space buffer (on device).
  real *rbuf;                // Real-space buffer (on device).
  real normFac;              // Normalization.
  bool forwardNorm;          // Normalize in r2c (true) or c2r (false) FFT.
  mugy_cufft_plan plan_r2c, plan_c2r;  // Plans.
  mugy_custream stream;
};

struct mugy_fft_fam_dev {
  struct mugy_fft_dev *xy;    // 2D FFT of x-y planes.
  struct mugy_fft_dev *xy_a;  // 2D FFT of (aliased) x-y planes.
};
