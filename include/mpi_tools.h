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

extern int myRank, totNumProcs;  // Rank of this process & total number of processes.
extern int numProcs[nDim+1];     // Processes along x,y,z and species.
extern MPI_Comm cartComm;        // Cartesian communicator.
extern int cartRank;             // Rank in cartCOMM.
extern MPI_Comm *sub1dComm;      // 1D subcommunicators along each direction.
extern int sub1dRank[nDim+1];    // ID (rank) in the 1D xpec,Z,X,Y subcommunicators.

void init_mpi(int argc, char *argv[]);  // Initialize MPI.

// Initialize sub-communicators.
void init_comms(struct grid grid, struct speciesParameters spec);

// Allocate array var and copy numE elements to it from src.
void allocAndCopyVar_int(int **var, int *src, const int numE);
void allocAndCopyVar_real(real **var, real *src, const int numE);

// Distribute degrees of freedom amongst MPI processes in 1D.
void distribute1dDOFs(const int procs, const int procID, const int globalDOFs, int *localDOFs, int *firstDOF);

// Distribute s,Z,X,Y degrees of freedom amongst MPI processes.
void distributeDOFs(struct grid globalGrid, struct speciesParameters globalSpec, struct grid *localGrid, struct speciesParameters *localSpec);

void terminate_mpi();  // Terminate MPI.

#endif
