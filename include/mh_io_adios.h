/* mugy: mh_io_adios.h

   Header file for ADIOS IO module.
*/

#ifndef MUGY_IO_ADIOS
#define MUGY_IO_ADIOS

#include "adios2_c.h"
#include "mh_userFLAGS.h"
#include "mh_data.h"
#include "mh_utilities.h"

extern adios2_adios *ad_ctx;  // ADIOS context used throughout our IO.

// Variable engine for outputting momk.
extern adios2_variable *ad_momk_var;
extern adios2_engine *ad_momk_eng;

// Variable engine for outputting realArray.
extern adios2_variable *ad_arr_var;
extern adios2_engine *ad_arr_eng;
// Variable engine for outputting fourierArray.
extern adios2_variable *ad_arrk_var;
extern adios2_engine *ad_arrk_eng;

#if USE_SINGLE_PRECISION
#define adios_real adios2_type_float
#define adios_fourier adios2_type_float_complex
#else
#define adios_real adios2_type_double
#define adios_fourier adios2_type_double_complex
#endif

#define adios_mint adios2_type_int32_t

struct ad_fhandle {
  char *fname;           // File name.
  adios2_io *io;         // IO component
  adios2_variable *var;  // Variable (we may define more later).
  adios2_engine *eng;    // Engine (responsible for reading/writing).
  bool isVarReal;        // Indicate if var is real (fourier otherwise).
};

struct mugy_ioManager {
  struct ioSetup setup;      // IO setup (input file instructions).
  adios2_adios *ctx;      // ADIOS context used throughout our IO.
  struct ad_fhandle *files;  // Pointers to adios file handles.
  mint numfiles;             // Number of files.
};

// Start the IO interface.
void init_io();

// Create files for IO.
void setup_files(struct grid globalGrid, struct grid localGrid, struct population globalPop, struct population localPop);

// Output Fourier-space array..
void write_fourierMoments(struct fourierArray arrkIn);
void write_fourierArray(struct fourierArray arrkIn);
void write_realArray(struct realArray arrIn);

// Finalize ADIOS IO.
void terminate_io();

#endif
