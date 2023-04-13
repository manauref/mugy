/* mugy: mh_io_adios.h
 *
 * Header file for ADIOS IO module.
 *
 */
#pragma once

#include "adios2_c.h"
#include "mh_userFLAGS.h"
#include "mh_comms.h"
#include "mh_data.h"
#include "mh_grid.h"
#include "mh_population.h"
#include "mh_utilities.h"

#ifdef USE_SINGLE_PRECISION
#define adios_real adios2_type_float
#define adios_fourier adios2_type_float_complex
#else
#define adios_real adios2_type_double
#define adios_fourier adios2_type_double_complex
#endif

#define adios_mint adios2_type_int32_t

// Container for IO instructions
struct mugy_ioSetup {
  char *inputFile;           // Name of input file.
  char *outputDir;           // Address of output directory.
  char *restartDir;          // Address of restart directory.
  bool isRestart;            // Is this simulation a restart of a previous one?
  bool outToOldDir;          // If restart, output to directory of previous run?
};

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
struct mugy_ioManager *mugy_io_init(struct mugy_comms comms); 

// Create a file holding global real(Fourier)Arrays.
struct mugy_ad_file *mugy_io_create_mugy_array_file(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid grid, enum mugy_datatype dtype);

// Create a file holding global real(Fourier) moments.
struct mugy_ad_file *mugy_io_create_moments_file(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid grid, struct mugy_population pop, enum mugy_datatype dtype);

// Create files for IO.
void mugy_io_setup_files(struct mugy_ioManager *ioman, struct mugy_grid grid, struct mugy_population pop);

// Output real(Fourier)-space array.
void mugy_io_write_mugy_array(struct mugy_ioManager *ioman, char* fname, struct mugy_ad_file *fh, struct mugy_array arr);

// Close a file given its mugy file handle.
void mugy_io_close_file(struct mugy_ad_file *fh);

// Finalize ADIOS IO.
void mugy_io_terminate(struct mugy_ioManager *ioman);
