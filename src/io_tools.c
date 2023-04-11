/* mugy: i0_tools.c
 *
 * Extra generic IO functions.
 *
 */

#include "mh_io_tools.h"
#include "mh_utilities.h"

void r0printf(char *str, mint rank) {
  // Print string out if this is the zeroth rank.
  if (rank == ioRank) printf(str);
}

// Print out elements in an array consecutively, separated by a space.
void arrPrintAll_mint(const mint *arr, mint numElements, const char *preStr, const char *postStr) {
  mint maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_mint(strOut, arr, numElements, preStr, postStr);
  printf(strOut);
}
void arrPrintAll_real(const real *arr, mint numElements, const char *preStr, const char *postStr) {
  mint maxCharacters = 9999;
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
