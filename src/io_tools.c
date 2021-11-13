/* mugy: i0_tools

   IO module.
*/

#include "io_tools.h"

void r0printf(char *str) {
  // Print string out if this is the zeroth rank.
  if (myRank == 0) printf(str);
}

// Print out elements in an array consecutively, separated by a space.
void arrPrint_int(const int *arr, const int numElements, const char *preStr, const char *postStr) {
  int maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_int(strOut, arr, numElements, preStr, postStr);
  r0printf(strOut);
}
void arrPrint_real(const real *arr, const int numElements, const char *preStr, const char *postStr) {
  int maxCharacters = 9999;
  char strOut[maxCharacters];
  arr2str_real(strOut, arr, numElements, preStr, postStr);
  r0printf(strOut);
}
