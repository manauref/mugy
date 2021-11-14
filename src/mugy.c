/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mugy.h"

int main(int argc, char *argv[]) {

  struct grid gridG, gridL;
  struct timeSetup timePars;
  struct speciesParameters specParsG, specParsL;
  struct ioSetup myIO;
  struct fieldParameters fieldPars; 

  struct realMoments mom, moma;
  struct fourierMoments momk, momka;

  init_mpi(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n" );

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &myIO, &gridG, &timePars, &specParsG, &fieldPars);
  init_comms(gridG, specParsG);

  // Set the number of cells in Fourier space and aliased real space.
  init_global_grids(&gridG);
  distributeDOFs(gridG, specParsG, &gridL, &specParsL);

//  printf(" Number of time steps and frames:           Nt       =%8d   |  nFrames  =%6d\n", 1000, timePars.nFrames);

  resource onResource = hostOnly;
  alloc_realMoments( gridL.fG.dual, specParsG, onResource, &mom);
  alloc_realMoments(gridL.fGa.dual, specParsG, onResource, &moma);
  alloc_fourierMoments( gridL.fG, specParsG, onResource, &momk);
  alloc_fourierMoments(gridL.fGa, specParsG, onResource, &momka);

  free_realMoments(&mom, onResource);
  free_realMoments(&moma, onResource);
  free_fourierMoments(&momk, onResource);
  free_fourierMoments(&momka, onResource);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.
  free_grid(&gridL);
  free_speciesPars(&specParsL);
  free_grid(&gridG);
  free_speciesPars(&specParsG);

  terminate_mpi();  // Finalize MPI.

  return 0;
}
