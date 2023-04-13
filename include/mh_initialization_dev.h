/* mugy: mh_initialization_dev.h
 *
 * Functions used to initialize devices (GPUs).
 *
 */
#pragma once

#include "mh_macros.h"
#include "mh_comms.h"

// Run the device initialization.
void device_init_dev(struct mugy_comms *comms);
