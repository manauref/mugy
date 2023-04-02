/* mugy: mh_userFLAGS.h

   Pre-processor flags for mugy
*/

#ifndef MUGY_PARAMETERS
#define MUGY_PARAMETERS

/* Specify whether to do a single or a double precision simulation.
     = 0 use double precision
     > 0 use single precision. */
#define USE_SINGLE_PRECISION 1 

// Time steppers supported. Options must be consistent with mh_macros.h
//   =4 Runge-Kutta 4th order.
#define TIME_STEPPER 4


#endif
