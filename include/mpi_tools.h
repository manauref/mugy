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

extern int myRank, numProcs;  // Rank of this process & total number of processes.
extern int xsProcs[4];        // Processes along x,y,z and species.
extern MPI_Comm cartCOMM;     // Cartesian communicator.

void init_mpi(int argc, char *argv[]);  // Initialize MPI.

// Initialize sub-communicators.
void init_comms(struct grid grid, struct speciesParameters spec);

void terminate_mpi();  // Terminate MPI.

#endif
