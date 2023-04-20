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
#include "mh_field.h"
#include "mh_utilities.h"
#include "mh_time.h"

#ifdef USE_SINGLE_PRECISION
#define MUGY_ADIOS_REAL adios2_type_float
#define MUGY_ADIOS_FOURIER adios2_type_float_complex
#else
#define MUGY_ADIOS_REAL adios2_type_double
#define MUGY_ADIOS_FOURIER adios2_type_double_complex
#endif

#define MUGY_ADIOS_MINT adios2_type_int32_t
#define MUGY_ADIOS_STRING adios2_type_string

// Container for IO instructions
struct mugy_io_pars {
  char *inputFile;           // Name of input file.
  char *outputDir;           // Address of output directory.
  char *restartDir;          // Address of restart directory.
  bool isRestart;            // Is this simulation a restart of a previous one?
  bool outToOldDir;          // If restart, output to directory of previous run?
};

struct mugy_ad_file {
  char *fname;              // File name.
  adios2_io *io;            // IO component
  adios2_variable **varG;   // Array of global variables i file.
  adios2_variable *time;    // Current simulation time.
  adios2_variable *frame;   // Frame number.
  adios2_engine *eng;       // Engine (responsible for reading/writing).
};

struct mugy_io {
  struct mugy_io_pars setup;   // IO setup (input file instructions).
  adios2_adios *ctx;      // ADIOS context used throughout our IO.
  struct mugy_ad_file **files; // Pointers to adios file handles.
  mint numfiles;          // Number of files.
  mint world_rank;        // Used so that only one rank writes something.
  // Buffers to write global single variables.
  double time_buf[1];
  mint frame_buf[1];
};

// Start the IO interface.
struct mugy_io *mugy_io_init(struct mugy_comms *comms); 

// Create a file holding global real(Fourier)Arrays.
struct mugy_ad_file *mugy_io_create_mugy_array_file(struct mugy_io *ioman, char* fname,
  struct mugy_grid *grid, enum mugy_data_types dtype, enum mugy_grid_types gtype);

// Create a file holding global real(Fourier) moments.
struct mugy_ad_file *mugy_io_create_moments_file(struct mugy_io *ioman, char* fname,
  struct mugy_grid *grid, struct mugy_population *pop, enum mugy_data_types dtype, enum mugy_grid_types gtype);

// Create a file for a mugy_array holding ncomp quantities per species on an perpendicular plane.
struct mugy_ad_file *mugy_io_create_population_perp_file(struct mugy_io *ioman, char* fname,
  struct mugy_grid *grid, struct mugy_population *pop,
  enum mugy_data_types dtype, enum mugy_grid_types gridtype, mint ncomp, mint zIdx);

// Create a file storing fourier-space moments and fields to restart the simulation.
struct mugy_ad_file *mugy_io_create_restart_file(struct mugy_io *ioman,
  struct mugy_grid *grid, struct mugy_population *pop);

// Create files for IO.
void mugy_io_setup_files(struct mugy_io *ioman, struct mugy_grid *grid, struct mugy_population *pop);

// Output real(Fourier)-space array.
void mugy_io_write_mugy_array(struct mugy_io *ioman, char* fname, struct mugy_ad_file *fh,
  struct mugy_array *arr, mint frame, double time);

// Write out a mugy array that evolves in time.
void mugy_io_write_mugy_array_step(struct mugy_io *ioman, char* fname, struct mugy_ad_file *fin,
  struct mugy_array *arr, mint frame, double time);

// Write out a frame of all dynamical quantities of interest.
void mugy_io_write_frame(struct mugy_io *ioman, struct mugy_population *pop,
  struct mugy_field *field, mint stepIdx, mint frame, double time);

// Write out the file used to restart a simulation.
void mugy_io_write_restart(struct mugy_io *ioman, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_time *time, mint stepIdx);

// Close a file given its mugy file handle.
void mugy_io_close_file(struct mugy_ad_file *fh);

// Finalize ADIOS IO.
void mugy_io_terminate(struct mugy_io *ioman);
