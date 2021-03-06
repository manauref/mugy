/* mugy: io_tools.h

   IO module header file.
*/

#ifndef IO_TOOLS 
#define IO_TOOLS 

#include <string.h>   // e.g. for memcpy.
#include "utilities.h"
#include "data_mugy.h"
#include "mpi_tools.h"
#include "io_adios.h"

// Print string out if this is the zeroth rank.
void r0printf(char *str);

// Print out elements in an array consecutively, separated by a space.
void arrPrint_mint(const mint *arr, const mint numElements, const char *preStr, const char *postStr);
void arrPrint_real(const real *arr, const mint numElements, const char *preStr, const char *postStr);

//void writeFields(); 

#endif
