/* A mugy header file:

   Utility functions used in mugy.
*/

#ifndef INITIALIZATION
#define INITIALIZATION

#include "parameters.h"
#include "data_mugy.h"
#include "utilities.h"
#include "alloc_mugy.h"
#include <string.h>   // e.g. for strcat, strlen.
#include "mpi_tools.h"
#include "io_tools.h"

// Read a real variable from input file. Need this function to support
// reading floats and doubles with the same function call.
void readFileVar_real(FILE *fp, const int numElements, real *var);
void readFileVar_int(FILE *fp, const int numElements, int *var);
void allocAndReadFileVar_real(FILE *fp, const int numElements, real **var);
void allocAndReadFileVar_int(FILE *fp, const int numElements, int **var);

// Read input values from input file.
void read_inputFile(const char *fileNameIn, struct grid *grid, struct timeSetup *time,
                    struct speciesParameters *spec, struct fieldParameters *field);

// Read command line arguments and input file.
void read_inputs(int argc, char *argv[], struct ioSetup *ioSet, struct grid *grid, struct timeSetup *time,
                 struct speciesParameters *spec, struct fieldParameters *field);

// Set number of cells in de-aliased, aliased and real space global grids.
void init_global_grids(struct grid *grid);

// Allocate various fields needed.
void allocate_fields(struct grid localGrid, struct speciesParameters localSpec);

// Impose the initial conditions on the moments and the electrostatic potential.
void set_initialCondition(struct grid localGrid, struct speciesParameters localSpec);

// Deallocate fields.
void free_fields();

void free_grid(struct grid *grid); // Free arrays in grids.
void free_speciesPars(struct speciesParameters *spec);  // Free arrays in species parameters.
#endif
