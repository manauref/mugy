/* mugy: mh_comms.h
 *
 * MPI-related operations.
 *
 */
#pragma once

#include "mh_macros.h"
#include "mh_grid.h"
#include "mh_population.h"

struct mugy_comms_sub {
  void* comm;       // MPI/NCCL communicator.
  mint dim, size;   // Dimensionality and size of the communicator.
  mint *decomp;     // Decomposition (number of processes) in each direction;
  mint rank;        // Rank ID within this communicator.
  mint *coord;      // Coordinates of this rank in the grid of processes.
};

struct mugy_comms {
  struct mugy_comms_sub world;
  struct mugy_comms_sub *sub1d, *sub2d, *sub3d, *sub4d; // 1D, 2D, 3D, 4d subcomms.
};

struct mugy_comms *mugy_comms_init(mint argc, char *argv[]);  // Initialize MPI.

// Initialize sub-communicators.
void mugy_comms_sub_init(struct mugy_comms *comms, struct mugy_grid *grid, struct mugy_population *pop);

// Distribute s,Z,X,Y degrees of freedom amongst MPI processes.
void mugy_comms_distributeDOFs(struct mugy_comms *comms, struct mugy_grid *grid, struct mugy_population *pop);

void mugy_comms_terminate(struct mugy_comms *comms);  // Terminate communications.
