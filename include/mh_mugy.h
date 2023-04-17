/* mugy: mh_mugy.h
 *
 * mugy (header file)
 * A gyrofluid code for multiscale simulations.

 * This file contains global macros and inclusion of various header files.
 *
 */

// Pre-processor flags.
#include "mh_userFLAGS.h"

// Generic utility functions (with few or no dependencies).
#include "mh_utilities.h"

// MPI-related infrastructure.
#include "mh_comms.h"

// IO module.
#include "mh_io_tools.h"

// Data structures specific to mugy.
#include "mh_data.h"

// Grid object and its methods.
#include "mh_grid.h"

// Population object and functions to manipulate it.
#include "mh_population.h"

// Field object and methods to query/modify it.
#include "mh_field.h"

// Functions that allocate mugy data.
#include "mh_alloc.h"

// FLR and Linear operators.
#include "mh_flr.h"

// Functions used to initialize the simulation.
#include "mh_initialization.h"

// FFT operators.
#include "mh_ffts.h"

// Time integration objects and functions.
#include "mh_time.h"
