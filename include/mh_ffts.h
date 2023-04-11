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
#pragma once

#include "mh_comms.h"
#include "mh_array.h"
#include "mh_grid.h"

typedef struct mugy_ffts mugy_ffts;

struct mugy_ffts *fft_init(struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_comms comms);

void fft_xy_c2r(struct mugy_ffts *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum resource_comp res);
void fft_xy_r2c(struct mugy_ffts *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum resource_comp res);

void fft_terminate(struct mugy_ffts *ffts);
