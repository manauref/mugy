/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mugy.h"

int main(int argc, char *argv[]) {

  int numSpecies = 1;
  int numMoments = 2;
  struct gridType grid;
  struct ioSetup myIO;

  struct realMoments mom, moma;
  struct fourierMoments momk, momka;

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &myIO, &grid);

  // Set the number of cells in Fourier space and aliased real space.
  init_gridsG(grid.NkxG, &grid);

  resource onResource = hostOnly;
  alloc_realMoments(numSpecies, numMoments, grid, onResource, &mom);
  alloc_fourierMoments(numSpecies, numMoments, grid, onResource, &momk);
  alloc_aliasedRealMoments(numSpecies, numMoments, grid, onResource, &moma);
  alloc_aliasedFourierMoments(numSpecies, numMoments, grid, onResource, &momka);

  free_realMoments(&mom, onResource);
  free_fourierMoments(&momk, onResource);
  free_realMoments(&moma, onResource);
  free_fourierMoments(&momka, onResource);

  return 0;
}
