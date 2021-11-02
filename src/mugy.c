/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include <string.h>   // e.g. for strcat.
#include "mugy.h"

int main(int argc, char *argv[]) {

  int numSpecies = 1;
  int numMoments = 2;
  struct gridType grid;
  char *inputFile;           // Name of input file.
  char *outputDir;           // Address of output directory.
  char *restartDir;          // Address of restart directory.
  bool isRestart;            // Is this simulation a restart of a previous one?
  bool outToOldDir;          // If restart, output to directory of previous run?
  char *checkFile;

  real *mom, *mom_a;
  fourier *momk, *momk_a;

  // Check for commandline arguments.
  // Currently we expect two arguments in this order:
  //   1) name of the input file.
  //   2) absolute address of the output directory.
  // For restarts add the restart directory as a 3rd argument.
  if (argc < 3) {
    printf("\n --> Not enough inputs. Need input file and output directory as command line arguments.\n");
    abortSimulation(" TERMINATING ");
  } else {
    inputFile = argv[1];
    outputDir = argv[2];
    isRestart   = false;
    outToOldDir = false;
    if (argc > 3) {  // Simulation is a restart of a previous one.
      restartDir = argv[3];
      isRestart  = true;
      // Output to the same directory as the previous run? 
      checkFile = malloc(strlen(outputDir)+strlen("phik.bp"));
      checkFile[0] = '\0';
      strcat(strcat(checkFile,outputDir),"phik.bp");
      fileExists(checkFile);
      free(checkFile);
    }
  }
  
  read_inputs(inputFile, grid.NkxG, grid.kxMin);

  // Set the number of cells in Fourier space and aliased real space.
  init_gridsG(grid.NkxG, grid);

  mom  = allocRealMoments(numSpecies, numMoments, grid);
  momk = allocFourierMoments(numSpecies, numMoments, grid);
  mom_a  = allocAliasedRealMoments(numSpecies, numMoments, grid);
  momk_a = allocAliasedFourierMoments(numSpecies, numMoments, grid);

  free(mom);
  free(momk);
  free(mom_a);
  free(momk_a);

  return 0;
}
