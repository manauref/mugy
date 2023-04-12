/* mugy: mh_userFLAGS.h
 *
 * Pre-processor flags for mugy
 * 
 */
#pragma once

/* Specify whether to do a single or a double precision simulation.
     #define: use single precision
     #undef:  use double precision. */
#undef USE_SINGLE_PRECISION

// Time steppers supported. Options must be consistent with mh_macros.h
//   =4 Runge-Kutta 4th order.
#define TIME_STEPPER 4


