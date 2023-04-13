/* mugy: mh_grid.h
 *
 * Grid object and related methods.
 *
 */
#pragma once

#include "mh_macros.h"

struct mugy_realGrid {
  mint Nx[nDim];        // Number of cells.
  mint NxTot;           // Total number of cells.
  mint NxyTot;          // Total number of cells in an x-y plane.
  real xMin[nDim];      // Minimum coordinates.
  real xMax[nDim];      // Maximum coordinates.
  real dx[nDim];        // Cell length.
  real *x;              // Coordinates in each direction.
  mint globalOff[nDim]; // Offset of first element in this process within the global domain.
};

struct mugy_fourierGrid {
  mint Nkx[nDim];        // Number of distinct absolute amplitude wavenumbers (counting k=0).
  mint Nekx[nDim];       // Number of elements in a k-space array (counting k=0 and negative k's).
  mint NekxTot;          // Total number of elements.
  mint NekxyTot;         // Total number of cells in an kx-ky plane.
  real kxMin[nDim];      // Minimum finite absolute amplitude wavenumbers.
  real *kx;              // Coordinates along each direction.
  struct mugy_realGrid dual;  // Real grid dual to this Fourier grid.
  real kxMaxDyn[nDim];   // Maximum k evolved. Above this we multiply time rates of change by zero.
  mint globalOff[nDim];  // Offset of first element in this process within the global domain.
};

struct mugy_grid_ada {
  struct mugy_fourierGrid al;    // Aliased Fourier grid.
  struct mugy_fourierGrid deal;  // Dealiased Fourier grid.
};

struct mugy_grid {
  struct mugy_grid_ada global;  // Global grids.
  struct mugy_grid_ada local;   // Local grids.
  mint mpiProcs[nDim];    // Number of MPI processes along each direction.
};

// Allocate mugy_grid object.
struct mugy_grid *mugy_grid_alloc();

// Initialize the global grids.
void mugy_grid_init_global(struct mugy_grid *grid, mint rank);

// Linear index given the nDim-dimensional subscript in a real/Fourier grid.
mint sub2lin_real(mint *xI, const struct mugy_realGrid grid);
mint sub2lin_fourier(mint *kxI, const struct mugy_fourierGrid grid);

// nDim-dimensional subscript given the linear index in a real/Fourier grid.
void lin2sub_real(mint *xI, mint lin, const struct mugy_realGrid grid);
void lin2sub_fourier(mint *kxI, mint lin, const struct mugy_fourierGrid grid);

// (x,y,z)/(kx,ky,kz) coordinates given the multidimensional xI/kxI index.
void get_x(real *x, mint *xI, const struct mugy_realGrid grid);
void get_kx(real *kx, mint *kxI, const struct mugy_fourierGrid grid);

// Deallocate memory used by grids.
void mugy_grid_free(struct mugy_grid *grid);
