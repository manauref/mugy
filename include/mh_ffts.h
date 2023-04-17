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

// Type of FFTs supported:
//   xy:       xy FFT of single 3D array on dealiased grids.
//   xy_a:     xy FFT of single 3D array on aliased grids.
//   mom_xy:   xy FFT of all moments on dealiased grids.
//   mom_xy_a: xy FFT of all moments on aliased grids.
enum mugy_fft_type {mugy_fft_xy, mugy_fft_mom_xy, mugy_fft_xy_a, mugy_fft_mom_xy_a};

typedef struct mugy_fft mugy_fft;

struct mugy_fft *mugy_fft_init(struct mugy_grid *grid, struct mugy_population_species *popL, struct mugy_comms *comms);

void mugy_fft_c2r(struct mugy_fft *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum mugy_fft_type ttype, enum mugy_resource_calc res);
void mugy_fft_r2c(struct mugy_fft *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum mugy_fft_type ttype, enum mugy_resource_calc res);

void mugy_fft_terminate(struct mugy_fft *ffts);
