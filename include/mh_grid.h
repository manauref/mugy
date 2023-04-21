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

// Linear index given the dimIn-dimensional subscript in a grid.
unsigned long mugy_grid_sub2lin(mint *xI, const struct mugy_grid_basic *grid, mint dimIn);

// Given the linear index 'lin', return the dimIn-dimensional index (subscript) xI.
void mugy_grid_lin2sub(mint *xI, unsigned long lin, const struct mugy_grid_basic *grid, mint dimIn);

// Obtain the coordinates given the dimIn-dimensional xI index.
void mugy_grid_get_x(real *x, mint *xI, const struct mugy_grid_basic *grid, mint dimIn);

// Deallocate memory used by grids.
void mugy_grid_free(struct mugy_grid *grid);
