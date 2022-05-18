/* mugy: io_adios.h

   Header file for ADIOS IO module.
*/

#ifndef IO_ADIOS
#define IO_ADIOS

#include "adios2_c.h"
#include "mh_userFLAGS.h"
#include "mh_data.h"
#include "mh_utilities.h"

extern adios2_adios *ad_ctx;  // ADIOS context used throughout our IO.

// Variable engine for outputting momk.
extern adios2_variable *ad_momk_var;
extern adios2_engine *ad_momk_eng;

#if USE_SINGLE_PRECISION > 0
#define adios_real adios2_type_float
#define adios_fourier adios2_type_float_complex
#else
#define adios_real adios2_type_double
#define adios_fourier adios2_type_double_complex
#endif

// Start the IO interface.
void init_io();

// Create files for IO.
void setup_files(struct grid globalGrid, struct grid localGrid, struct population globalPop, struct population localPop);

// Output Fourier-space moments.
void writeMoments_fourier(struct fourierMoments momkIn);

// Finalize ADIOS IO.
void terminate_io();

#endif
