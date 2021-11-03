/* 
   mugy (header file)
   A reduced gyrofluid code for ITG-ETG multiscale simulations.

   This file contains global macros and inclusion of various header files.
*/

/* .................. GLOBAL MACROS ..................... */

/* ............... END OF GLOBAL MACROS ................. */

// Pre-processor flags.
#include "mugyFLAGS.h"

// Global parameters.
#include "parameters.h"

// Generic utility functions (with few or no dependencies).
#include "utilities.h"

// Data structures specific to mugy.
#include "data_mugy.h"

// Functions that allocate mugy data.
#include "alloc_mugy.h"

// Functions used to initialize the simulation.
#include "initialization.h"
