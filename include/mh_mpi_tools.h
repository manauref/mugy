/* mugy: mh_mpi_tools.h

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

extern mint myRank, totNumProcs;  // Rank of this process & total number of processes.
extern mint numProcs[nDim+1];     // Processes along x,y,z and species.
extern MPI_Comm cartComm;         // Cartesian communicator.
extern mint cartRank;             // Rank in cartCOMM.

extern MPI_Comm *sub1dComm;     // 1D subcommunicators along each direction.
extern mint sub1dRank[nDim+1];  // ID (rank) in the 1D spec,Z,X,Y subcommunicators.
extern MPI_Comm *sComm, *zComm, *xComm, *yComm;  // Pointers for 1D comms.
extern mint sRank, zRank, xRank, yRank;          // Pointers for 1D rank IDs.

extern MPI_Comm *sub2dComm;   // 2D subcommunicators (e.g. XY comm).
extern mint sub2dRank[nDim];  // ID (rank) in the 2D (xy,xz,yz) subcommunicators.
extern MPI_Comm *xyComm, *xzComm, *yzComm;  // Pointers for 2D comms.
extern mint xyRank, xzRank, yzRank;         // Pointers for 2D rank IDs.

void init_mpi(mint argc, char *argv[]);  // Initialize MPI.

// Initialize sub-communicators.
void init_comms(struct mugy_grid grid, struct mugy_population pop);

// Distribute degrees of freedom amongst MPI processes in 1D.
void distribute1dDOFs(const mint procs, const mint procID, const mint globalDOFs, mint *localDOFs, mint *firstDOF);

// Distribute s,Z,X,Y degrees of freedom amongst MPI processes.
void distributeDOFs(struct mugy_grid globalGrid, struct mugy_population globalPop, struct mugy_grid *localGrid, struct mugy_population *localPop);

void terminate_mpi();  // Terminate MPI.

#endif
