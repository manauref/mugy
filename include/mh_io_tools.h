/* mugy: mh_io_tools.h

   IO module header file.
*/

#ifndef MUGY_IO_TOOLS 
#define MUGY_IO_TOOLS 

#include <string.h>   // e.g. for memcpy.
#include "mh_utilities.h"
#include "mh_data.h"
#include "mh_mpi_tools.h"
#include "mh_io_adios.h"

// Print string out if this is the zeroth rank.
void r0printf(char *str);

// Print out elements in an array consecutively, separated by a space.
void arrPrint_mint(const mint *arr, const mint numElements, const char *preStr, const char *postStr);
void arrPrint_real(const real *arr, const mint numElements, const char *preStr, const char *postStr);

//void writeFields(); 

#endif
