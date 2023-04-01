/* mugy: mh_mpi_data.h

   MPI-related data/parameters.
*/

#ifndef MUGY_MPI_DATA
#define MUGY_MPI_DATA

extern mint myRank, numProcs;  // Rank of this process & total number of processes.
extern mint xsProcs[4];        // Processes along x,y,z and species.
extern MPI_Comm cartCOMM;     // Cartesian communicator.

#endif
