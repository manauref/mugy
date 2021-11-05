/* mugy: mpi_tools

   MPI-related operations.
*/

#include "mpi_tools.h"

int myRank, numProcs;  // Rank of this process & total number of processes.
int xsProcs[4];        // Processes along x,y,z and species.
MPI_Comm cartCOMM;     // Cartesian communicator.

void init_mpi(int argc, char *argv[]) {
  // Initialize MPI, get rank of this process and total number of processes.
  MPI_Init(&argc, &argv);  // Initialize MPI.
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  // Rank of this MPI process.
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);  // Number of MPI processes.
}

void init_comms(struct grid grid, struct speciesParameters spec) {
  // Initialize the various sub-communicators needed.
  
  // Check the number of MPI processes is correct.
  if (prod_int(grid.mpiProcs,3)*spec.mpiProcs != numProcs) {
    printf(" Number of MPI processes in input file (%d %d) differs from that in mpirun (%d).\n",
           prod_int(grid.mpiProcs,3),spec.mpiProcs, numProcs);
    //abortSimulation(" Terminating...\n");
  }

  memcpy(xsProcs, &grid.mpiProcs, 3*sizeof(int));
  xsProcs[3] = spec.mpiProcs;

  // MPI_CART boundary conditions. True=periodic.
  int cartCommBCs[4] = {true, true, true, false};
  // Let MPI assign arbitrary ranks.
  int reorder = true;

  // Create a 4D Cartesian communicator (3D space + species).
//  MPI_Cart_create(MPI_COMM_WORLD, 4, xsProcs, cartCommBCs, reorder, &cartCOMM);
}

void terminate_mpi() {
  // Close MPI up, thoroughly.
  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

//  MPI_Comm_free(&cartCOMM);

  MPI_Finalize();
}
