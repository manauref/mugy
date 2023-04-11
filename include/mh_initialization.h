/* mugy: mh_initialization.h
 *
 * Utility functions used in mugy.
 *
 */
#pragma once

#include "mh_data.h"
#include "mh_grid.h"
#include "mh_population.h"
#include "mh_field.h"
#include "mh_io_adios.h"
#include "mh_ffts.h"
#include "mh_timestepping.h"

// Read a real variable from input file. Need this function to support
// reading floats and doubles with the same function call.
void readFileVar_real(FILE *fp, const mint numElements, real *var);
void readFileVar_mint(FILE *fp, const mint numElements, mint *var);
// Read species parameter composed of numElements[s] for the s-th species.
void readFileSpeciesPar_mint(mint **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements);
void readFileSpeciesPar_real(real **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements);

void read_inputs(mint argc, char *argv[], struct mugy_ioSetup *ioSet, struct mugy_grid *grid, struct mugy_timeSetup *time,
                 struct mugy_population *pop, struct mugy_fieldParameters *field, mint rank);

void device_init(struct mugy_comms *comms);

// Allocate time dependent fields needed.
void allocate_dynfields(struct mugy_grid localGrid, struct mugy_population *localPop);

void init_global_grids(struct mugy_grid *globalGrid, mint rank);

void set_initialCondition(struct mugy_grid globalGrid, struct mugy_grid localGrid, struct mugy_population globalPop,
  struct mugy_population *localPop, struct mugy_ffts *fftMan, struct mugy_ioManager *ioman);

// Deallocate fields.
void free_fields();

void free_grid(struct mugy_grid *grid); // Free arrays in grids.
void free_population(struct mugy_population *pop);  // Free arrays in population struct.
