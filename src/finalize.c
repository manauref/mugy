/* mugy: finalize.c
   
   Function to end the simulation in a tidy manner.
*/
#include "finalize.h"

void terminate_all() {
  // Call the termination of the various mugy components.
  
  terminate_io();  // Close IO interface.

  terminate_mpi();  // Finalize MPI.

}
