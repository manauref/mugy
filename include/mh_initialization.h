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
                 struct mugy_population *pop, struct mugy_field *field, mint rank);

void device_init(struct mugy_comms *comms);

void set_initialConditions(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid,
  struct mugy_ffts *fftMan, struct mugy_ioManager *ioman);
