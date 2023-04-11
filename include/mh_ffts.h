/* mugy: mh_ffts.h
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
 */

#ifndef MUGY_FFTS
#define MUGY_FFTS

#include "mh_fftw_wrap.h"
#include "mh_comms.h"
#include "mh_array.h"
#include "mh_grid.h"
#include <stdbool.h>

// Info needed for a single FFT on the host.
struct mugy_fft_ho {
  mugy_fftw_fourier *kbuf;  // Fourier-space buffer.
  real *rbuf;               // Real-space buffer.
  real normFac;             // Normalization.
  bool forwardNorm;         // Normalize in r2c (true) or c2r (false) FFT.
  mugy_fftw_plan plan_r2c, plan_c2r;  // Plans.
};

typedef struct mugy_fft_dev mugy_fft_dev;

struct mugy_fft_hodev {
  struct mugy_fft_ho *ho;
  struct mugy_fft_dev *dev;
};

// Info needed for all FFTs in mugy.
struct mugy_ffts {
  struct mugy_fft_hodev xy;    // 2D FFT of x-y planes.
  struct mugy_fft_hodev xy_a;  // 2D FFT of (aliased) x-y planes.
};

void fft_init(struct mugy_ffts *ffts, struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_comms comms);

void fft_xy_c2r(struct mugy_ffts *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum resource_comp res);
void fft_xy_r2c(struct mugy_ffts *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum resource_comp res);

void fft_terminate(struct mugy_ffts *ffts);

#endif
