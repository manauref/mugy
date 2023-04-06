/* mugy: mh_io_adios.h

   Header file for ADIOS IO module.
*/

#ifndef MUGY_IO_ADIOS
#define MUGY_IO_ADIOS

#include "adios2_c.h"
#include "mh_userFLAGS.h"
#include "mh_data.h"
#include "mh_utilities.h"

#if USE_SINGLE_PRECISION
#define adios_real adios2_type_float
#define adios_fourier adios2_type_float_complex
#else
#define adios_real adios2_type_double
#define adios_fourier adios2_type_double_complex
#endif

#define adios_mint adios2_type_int32_t

struct mugy_ad_file {
  char *fname;           // File name.
  adios2_io *io;         // IO component
  adios2_variable *var;  // Variable (we may define more later).
  adios2_engine *eng;    // Engine (responsible for reading/writing).
  bool isVarReal;        // Indicate if var is real (fourier otherwise).
};

struct mugy_ioManager {
  struct mugy_ioSetup setup;   // IO setup (input file instructions).
  adios2_adios *ctx;      // ADIOS context used throughout our IO.
  struct mugy_ad_file **files; // Pointers to adios file handles.
  mint numfiles;          // Number of files.
};

// Start the IO interface.
void init_io(struct mugy_ioManager *ioman);

// Create a file holding global real(Fourier)Arrays.
struct mugy_ad_file *ad_create_file_realArray(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid);
struct mugy_ad_file *ad_create_file_fourierArray(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid);

// Create a file holding global real(Fourier) moments.
struct mugy_ad_file *ad_create_file_realMoments(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid, struct mugy_population globalPop, struct mugy_population localPop);
struct mugy_ad_file *ad_create_file_fourierMoments(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid, struct mugy_population globalPop, struct mugy_population localPop);

// Create files for IO.
void setup_files(struct mugy_ioManager *ioman, struct mugy_grid globalGrid, struct mugy_grid localGrid,
  struct mugy_population globalPop, struct mugy_population localPop);

// Output real(Fourier)-space array.
void write_realArray(struct mugy_ioManager *ioman, char* fname, struct mugy_ad_file *fhIn, struct mugy_realArray arrIn);
void write_fourierArray(struct mugy_ioManager *ioman, char* fname, struct mugy_ad_file *fhIn, struct mugy_fourierArray arrIn);

// Close a file given its mugy file handle.
void io_close_file(struct mugy_ad_file *fh);

// Finalize ADIOS IO.
void io_terminate(struct mugy_ioManager *ioman);

#endif
