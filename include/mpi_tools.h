/* mugy: mpi_tools header file.

   MPI-related operations.
*/

#ifndef MPI_TOOLS 
#define MPI_TOOLS 

#include <mpi.h>
#include <stdbool.h>  // e.g. for bool, true, false.
#include <string.h>   // e.g. for memcpy.
#include "utilities.h"
#include "data_mugy.h"
#include "alloc_mugy.h"
#include "io_tools.h"

extern mint myRank, totNumProcs;  // Rank of this process & total number of processes.
extern mint numProcs[nDim+1];     // Processes along x,y,z and species.
extern MPI_Comm cartComm;        // Cartesian communicator.
extern mint cartRank;             // Rank in cartCOMM.
extern MPI_Comm *sub1dComm;      // 1D subcommunicators along each direction.
extern mint sub1dRank[nDim+1];    // ID (rank) in the 1D xpec,Z,X,Y subcommunicators.

void init_mpi(mint argc, char *argv[]);  // Initialize MPI.

// Initialize sub-communicators.
void init_comms(struct grid grid, struct population pop);

// Distribute degrees of freedom amongst MPI processes in 1D.
void distribute1dDOFs(const mint procs, const mint procID, const mint globalDOFs, mint *localDOFs, mint *firstDOF);

// Distribute s,Z,X,Y degrees of freedom amongst MPI processes.
void distributeDOFs(struct grid globalGrid, struct population globalPop, struct grid *localGrid, struct population *localPop);

void terminate_mpi();  // Terminate MPI.

#endif
