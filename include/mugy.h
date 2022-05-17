/* 
   mugy (header file)
   A reduced gyrofluid code for ITG-ETG multiscale simulations.

   This file contains global macros and inclusion of various header files.
*/

// Pre-processor flags.
#include "mugyFLAGS.h"

// Global parameters.
#include "parameters.h"

// Generic utility functions (with few or no dependencies).
#include "utilities.h"

// MPI-related infrastructure.
#include "mpi_tools.h"

// IO module.
#include "io_tools.h"

// Data structures specific to mugy.
#include "data_mugy.h"

// Functions that allocate mugy data.
#include "alloc_mugy.h"

// Functions used to initialize the simulation.
#include "initialization.h"

// Module in charge of closing mugy sim.
#include "finalize.h"
