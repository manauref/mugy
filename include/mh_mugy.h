/* 
   mugy (header file)
   A reduced gyrofluid code for ITG-ETG multiscale simulations.

   This file contains global macros and inclusion of various header files.
*/

// Pre-processor flags.
#include "mh_userFLAGS.h"

// Global parameters.
#include "mh_parameters.h"

// Generic utility functions (with few or no dependencies).
#include "mh_utilities.h"

// MPI-related infrastructure.
#include "mh_mpi_tools.h"

// IO module.
#include "mh_io_tools.h"

// Data structures specific to mugy.
#include "mh_data.h"

// Functions that allocate mugy data.
#include "mh_alloc.h"

// Functions used to initialize the simulation.
#include "mh_initialization.h"

// Module in charge of closing mugy sim.
#include "mh_finalize.h"
