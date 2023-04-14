/* mugy: field.c
 *
 * Field object and operations to compute it.
 *
 */
#include "mh_field.h"
#include <stdlib.h>  // for malloc.

struct mugy_field *mugy_field_alloc() {
  // Allocate the field object.
  struct mugy_field *field = (struct mugy_field *) malloc(sizeof(struct mugy_field));
  return field;
}

void mugy_field_init(struct mugy_field *field, struct mugy_grid *grid, struct mugy_population *pop) {
#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Allocate Fourier-space potential.
  field->phik = mugy_array_alloc(fourier_enum, grid->local.deal.NekxTot, onResource);

  // Allocate Fourier-space gyroaveraged potentials, 3 for each species:
  //   <phi>, 0.5*hatLap <phi> and 0.5*hathatLap <phi>.
  field->gyrophik = mugy_array_alloc(fourier_enum, 3 * pop->local.numSpecies * grid->local.deal.NekxTot, onResource);
}

void mugy_field_free(struct mugy_field *field) {
#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Free potentials.
  mugy_array_free(field->phik, onResource); 
  mugy_array_free(field->gyrophik, onResource); 

  // Free the field object.
  free(field);
}
