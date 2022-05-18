/* mugy: mpi_data.h

   MPI-related data/parameters.
*/

#ifndef MPI_DATA
#define MPI_DATA

extern mint myRank, numProcs;  // Rank of this process & total number of processes.
extern mint xsProcs[4];        // Processes along x,y,z and species.
extern MPI_Comm cartCOMM;     // Cartesian communicator.

#endif
