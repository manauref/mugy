/* mugy: mh_comms.h

   MPI-related operations.
*/

#ifndef MUGY_MPI_TOOLS 
#define MUGY_MPI_TOOLS 

#include <mpi.h>
#include <stdbool.h>  // e.g. for bool, true, false.
#include "mh_utilities.h"
#include "mh_data.h"
#include "mh_alloc.h"
#include "mh_io_tools.h"

struct mugy_comms_sub {
  MPI_Comm comm;    // MPI communicator.
  mint dim, size;   // Dimensionality and size of the communicator.
  mint *decomp;     // Decomposition (number of processes) in each direction;
  mint rank;        // Rank ID within this communicator.
  mint *coord;      // Coordinates of this rank in the grid of processes.
};

struct mugy_comms {
  struct mugy_comms_sub world;
  struct mugy_comms_sub *sub1d, *sub2d, *sub3d, *sub4d; // 1D, 2D, 3D, 4d subcomms.
};

void comms_init(struct mugy_comms *comms, mint argc, char *argv[]);  // Initialize MPI.

// Initialize sub-communicators.
void comms_sub_init(struct mugy_comms *comms, struct mugy_grid grid, struct mugy_population pop);

// Distribute degrees of freedom amongst MPI processes in 1D.
void distribute1dDOFs(const mint procs, const mint procID, const mint globalDOFs, mint *localDOFs, mint *firstDOF);

// Distribute s,Z,X,Y degrees of freedom amongst MPI processes.
void distributeDOFs(struct mugy_comms comms, struct mugy_grid globalGrid, struct mugy_population globalPop, struct mugy_grid *localGrid, struct mugy_population *localPop);

void comms_terminate(struct mugy_comms *comms);  // Terminate communications.

#endif