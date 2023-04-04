/* mugy: finalize.c
   
   Function to end the simulation in a tidy manner.
*/
#include "mh_finalize.h"
#include "mh_ffts.h"

void terminate_all() {
  // Call the termination of the various mugy components.

  terminate_ffts();
  
  terminate_mpi();  // Finalize MPI.

}
