/* mugy: mh_grid.h
 *
 * Grid object and related methods.
 *
 */
#pragma once

#include "mh_macros.h"

// Types of grid used in mugy. We use the preprocessor in order
// to consistently generate strings for each type (see IO code).
#define MUGY_FOREACH_GRIDTYPE(GRIDTYPE) \
        GRIDTYPE(MUGY_REAL_GRID)   \
        GRIDTYPE(MUGY_FOURIER_GRID)  \

#define MUGY_GENERATE_ENUM(ENUM) ENUM,
enum mugy_grid_types {
    MUGY_FOREACH_GRIDTYPE(MUGY_GENERATE_ENUM)
};

struct mugy_grid_basic {
  mint Nx[nDim];        // Number of cells/coordinates.
  mint NxNonNeg[nDim];  // Number of nonnegative coordinates.
  mint NxTot;           // Total number of cells.
  mint NxyTot;          // Total number of cells in an x-y plane.
  real xMin[nDim];      // Minimum coordinates.
  real xMax[nDim];      // Maximum coordinates.
  real dx[nDim];        // Cell length.
  real *x;              // Coordinates in each direction.
  real *xperpSq;        // x^2+y^2.
  mint globalOff[nDim]; // Offset of first element in this process within the global domain.
  enum mugy_grid_types type;  // Type of the grid (real/fourier).
};

struct mugy_grid_chart {
  struct mugy_grid_basic *real;       // Fourier grid (dealiased).
  struct mugy_grid_basic *fourier;    // Fourier grid (dealiased).
  struct mugy_grid_basic *realAl;     // Aliased Fourier grid.
  struct mugy_grid_basic *fourierAl;  // Aliased Fourier grid.
};

struct mugy_grid {
  struct mugy_grid_chart *global;  // Global grids.
  struct mugy_grid_chart *local;   // Local grids.
  mint mpiProcs[nDim];    // Number of MPI processes along each direction.
};

// Allocate mugy_grid object.
struct mugy_grid *mugy_grid_alloc();

// Initialize the global grids.
void mugy_grid_init_global(struct mugy_grid *grid, mint rank);

// Linear index given the nDim-dimensional subscript in a real/Fourier grid.
mint mugy_grid_sub2lin_real(mint *xI, const struct mugy_grid_basic *grid);
mint mugy_grid_sub2lin_fourier(mint *kxI, const struct mugy_grid_basic *grid);

// nDim-dimensional subscript given the linear index in a real/Fourier grid.
void mugy_grid_lin2sub_real(mint *xI, mint lin, const struct mugy_grid_basic *grid);
void mugy_grid_lin2sub_fourier(mint *kxI, mint lin, const struct mugy_grid_basic *grid);

// Like mugy_grid_lin2sub_fourier but on a perpendicular plane only.
void mugy_grid_lin2sub_fourier_perp(mint *kxI, mint lin, const struct mugy_grid_basic *grid);

// (x,y,z)/(kx,ky,kz) coordinates given the multidimensional xI/kxI index.
void mugy_grid_get_x(real *x, mint *xI, const struct mugy_grid_basic *grid);
void mugy_grid_get_kx(real *kx, mint *kxI, const struct mugy_grid_basic *grid);

// Obtain the kx=(kx,ky) coordinates given the multidimensional
// kxI index. Assume the flat grid.kx array is organized as {kx,ky,kz}.
void mugy_grid_get_kx_perp(real *kx, mint *kxI, const struct mugy_grid_basic *grid);

// Deallocate memory used by grids.
void mugy_grid_free(struct mugy_grid *grid);
