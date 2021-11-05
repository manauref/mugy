/* mugy: io_tools header file.

   IO module.
*/

#ifndef IO_TOOLS 
#define IO_TOOLS 

#include <string.h>   // e.g. for memcpy.
#include "utilities.h"
#include "data_mugy.h"
#include "mpi_tools.h"

// Print string out if this is the zeroth rank.
void r0printf(char *str);

// Print out elements in an array consecutively, separated by a space.
void arrPrint_int(const int *arr, const int numElements, const char *preStr, const char *postStr);
void arrPrint_real(const real *arr, const int numElements, const char *preStr, const char *postStr);

#endif
