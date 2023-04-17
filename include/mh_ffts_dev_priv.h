/* mugy: mh_ffts_dev_priv.h
 * 
 * Private declarations for device FFTs.
 *
 */
#pragma once

#include "mh_macros.h"
#include <cufftXt.h>
#include "mh_cufft_utils.h"

#ifdef USE_SINGLE_PRECISION
typedef cufftComplex mugy_cufft_fourier_t;

const cudaDataType MUGY_CUFFT_REAL = CUDA_R_32F;
const cudaDataType MUGY_CUFFT_FOURIER = CUDA_C_32F;
const cudaDataType MUGY_CUFFT_EXEC_FOURIER = CUDA_C_32F;
#else
typedef cufftDoubleComplex mugy_cufft_fourier_t;

const cudaDataType MUGY_CUFFT_REAL = CUDA_R_64F;
const cudaDataType MUGY_CUFFT_FOURIER = CUDA_C_64F;
const cudaDataType MUGY_CUFFT_EXEC_FOURIER = CUDA_C_64F;
#endif

typedef cufftHandle mugy_cufft_plan;
typedef cudaStream_t mugy_custream;

// Info needed for a single FFT on the device.
struct mugy_fft_dev {
  mugy_cufft_fourier_t *kbuf;  // Fourier-space buffer (on device).
  real *rbuf;                // Real-space buffer (on device).
  real normFac;              // Normalization.
  bool forwardNorm;          // Normalize in r2c (true) or c2r (false) FFT.
  mugy_cufft_plan plan_r2c, plan_c2r;  // Plans.
  mugy_custream stream;
};

struct mugy_fft_fam_dev {
  // FFTs on dealised grids.
  struct mugy_fft_dev *xy;      // xy FFT of single 3D array.
  struct mugy_fft_dev *mom_xy;  // xy FFT of all moments.
  // FFTs on aliased grids.
  struct mugy_fft_dev *xy_a;      // xy FFT of single 3D array.
  struct mugy_fft_dev *mom_xy_a;  // xy FFT of all moments.
};
