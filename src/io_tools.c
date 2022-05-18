/* mugy: i0_tools

   IO module.
*/

#include "mh_io_tools.h"

void r0printf(char *str) {
  // Print string out if this is the zeroth rank.
  if (myRank == 0) printf(str);
}

// Print out elements in an array consecutively, separated by a space.
void arrPrint_mint(const mint *arr, const mint numElements, const char *preStr, const char *postStr) {
  mint maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_mint(strOut, arr, numElements, preStr, postStr);
  r0printf(strOut);
}
void arrPrint_real(const real *arr, const mint numElements, const char *preStr, const char *postStr) {
  mint maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_real(strOut, arr, numElements, preStr, postStr);
  r0printf(strOut);
}
