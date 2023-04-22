/* mugy: mh_io_utilities.h
 *
 * IO module header file.
 *
 */
#pragma once

#include "mh_macros.h"

// Print string out if this is the zeroth rank.
void r0printf(char *str, mint rank);

// Print out elements in an array consecutively, separated by a space.
void arrPrint_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr);
void arrPrint_real(const real *arr, mint numElements, const char *preStr, const char *postStr);

// Only let a single rank print.
void arrPrintS_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr, mint rank);
void arrPrintS_real(const real *arr, mint numElements, const char *preStr, const char *postStr, mint rank);

// Print a string followed by a single value and another string.
void valPrintS_mint(mint val, const char *preStr, const char *postStr, mint rank);
void valPrintS_real(real val, const char *preStr, const char *postStr, mint rank);
void valPrintS_char(char *val, const char *preStr, const char *postStr, mint rank);