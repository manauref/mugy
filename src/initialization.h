/* A mugy header file:

   Utility functions used in mugy.
*/

#ifndef INITIALIZATION
#define INITIALIZATION

#include "data_mugy.h"
#include "utilities.h"
#include <string.h>   // e.g. for strcat, strlen.

// Read input values from input file.
void read_inputFile(const char *fileNameIn, int *NkxG, real *kxMin);

// Read command line arguments and input file.
void read_inputs(int argc, char *argv[], struct ioSetup *ioSet, struct gridType *grid);

// Set number of cells in de-aliased, aliased and real space global grids.
void init_gridsG(const int *userNkxG, struct gridType *grid);

#endif
