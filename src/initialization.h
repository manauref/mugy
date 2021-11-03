/* A mugy header file:

   Utility functions used in mugy.
*/

#ifndef INITIALIZATION
#define INITIALIZATION

#include "data_mugy.h"
#include "utilities.h"

// Read input values from input file.
void read_inputs(const char *fileNameIn, int *NkxG, real *kxMin);

// Set number of cells in de-aliased, aliased and real space global grids.
void init_gridsG(const int *userNkxG, struct gridType *grid);

#endif
