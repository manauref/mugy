/* mugy: mh_io_tools.h
 *
 * IO module header file.
 *
 */

#ifndef MUGY_IO_TOOLS 
#define MUGY_IO_TOOLS 

#include <string.h>   // e.g. for memcpy.
#include "mh_utilities.h"
#include "mh_data.h"
#include "mh_comms.h"
#include "mh_io_adios.h"

// Print string out if this is the zeroth rank.
void r0printf(char *str, mint rank);

// Print out elements in an array consecutively, separated by a space.
void arrPrint_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr);
void arrPrint_real(const real *arr, mint numElements, const char *preStr, const char *postStr);

// Only let a single rank print.
void arrPrintS_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr, mint rank);
void arrPrintS_real(const real *arr, mint numElements, const char *preStr, const char *postStr, mint rank);

//void writeFields(); 

#endif
