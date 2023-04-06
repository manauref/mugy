/* mugy: mh_initialization.h

   Utility functions used in mugy.
*/

#ifndef MUGY_INITIALIZATION
#define MUGY_INITIALIZATION

#include "mh_data.h"

// Read a real variable from input file. Need this function to support
// reading floats and doubles with the same function call.
void readFileVar_real(FILE *fp, const mint numElements, real *var);
void readFileVar_mint(FILE *fp, const mint numElements, mint *var);
// Read species parameter composed of numElements[s] for the s-th species.
void readFileSpeciesPar_mint(mint **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements);
void readFileSpeciesPar_real(real **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements);

// Read input values from input file.
void read_inputFile(const char *fileNameIn, struct grid *grid, struct timeSetup *time,
                    struct population *pop, struct fieldParameters *field);

// Read command line arguments and input file.
void read_inputs(mint argc, char *argv[], struct ioSetup *ioSet, struct grid *grid, struct timeSetup *time,
                 struct population *pop, struct fieldParameters *field);

// Set number of cells in de-aliased, aliased and real space global grids.
void init_global_grids(struct grid *grid);

// Allocate time dependent fields needed.
void allocate_dynfields(struct grid localGrid, struct population *localPop);

// Impose the initial conditions on the moments and the electrostatic potential.
void set_initialCondition(struct grid localGrid, struct population *localPop, struct mugy_ioManager *ioman);

// Run the full initialization.
void init_all(mint argc, char *argv[], struct mugy_ioManager *ioman, struct grid *gridG, struct grid *gridL, struct timeSetup *timePars,
              struct population *popG, struct population *popL, struct fieldParameters *fieldPars);

// Deallocate fields.
void free_fields();

void free_grid(struct grid *grid); // Free arrays in grids.
void free_population(struct population *pop);  // Free arrays in population struct.
#endif
