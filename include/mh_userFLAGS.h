/* mugy: mh_userFLAGS.h

   Pre-processor flags for mugy
*/

#ifndef MUGY_PARAMETERS
#define MUGY_PARAMETERS

/* Specify whether to do a single or a double precision simulation.
     = false: use double precision
     = true: use single precision. */
#define USE_SINGLE_PRECISION true 

// Time steppers supported. Options must be consistent with mh_macros.h
//   =4 Runge-Kutta 4th order.
#define TIME_STEPPER 4


#endif
