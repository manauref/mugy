/* mugy: i0_tools.c
 *
 * Extra generic IO functions.
 *
 */

#include "mh_io_utilities.h"
#include "mh_utilities.h"
#include <string.h>

void r0printf(char *str, mint rank) {
  // Print string out if this is the zeroth rank.
  if (rank == ioRank) printf(str);
}

// Print out elements in an array consecutively, separated by a space.
void arrPrintAll_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr) {
  const mint maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_mint(strOut, arr, numElements, preStr, postStr);
  printf(strOut);
}
void arrPrintAll_real(const real *arr, mint numElements, const char *preStr, const char *postStr) {
  const mint maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_real(strOut, arr, numElements, preStr, postStr);
  printf(strOut);
}

// Only let a single rank print.
void arrPrintS_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr, mint rank) {
  if (rank == ioRank) arrPrintAll_mint(arr, numElements, preStr, postStr);
}
void arrPrintS_real(const real *arr, mint numElements, const char *preStr, const char *postStr, mint rank) {
  if (rank == ioRank) arrPrintAll_real(arr, numElements, preStr, postStr);
}

void valPrintS_mint(mint val, const char *preStr, const char *postStr, mint rank) {
  // Print a string followed by a single mint value and another string.
  const mint arr[] = {val};
  if (rank == ioRank) arrPrintAll_mint(arr, 1, preStr, postStr);
}
void valPrintS_real(real val, const char *preStr, const char *postStr, mint rank) {
  // Print a string followed by a single real value and another string.
  const real arr[] = {val};
  if (rank == ioRank) arrPrintAll_real(arr, 1, preStr, postStr);
}
void valPrintS_char(char *val, const char *preStr, const char *postStr, mint rank) {
  // Print a string followed by a single real value and another string.
  if (rank == ioRank) {
    char strOut[strlen(val)+strlen(preStr)+strlen(postStr)+1];
    strcpy(strOut, preStr);
    strcat(strcat(strOut,val),postStr);
    printf(strOut);
  }
}
