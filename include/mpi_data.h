/* mugy: mpi_data header file.

   MPI-related data/parameters.
*/

#ifndef MPI_DATA
#define MPI_DATA

extern int myRank, numProcs;  // Rank of this process & total number of processes.
extern int xsProcs[4];        // Processes along x,y,z and species.
extern MPI_Comm cartCOMM;     // Cartesian communicator.

#endif
