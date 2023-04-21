/* mugy
   
   A collection of utility functions.
*/
#include "mh_utilities.h"

mint prod_mint(const mint *arrIn, const mint numElements) {
  /* Compute the product of the elements in an array. */
  mint pOut = 1;
  for (mint i=0; i<numElements; i++) {pOut *= arrIn[i];}
  return pOut;
}
mint sum_mint(const mint *arrIn, const mint numElements) {
  /* Compute the sum of the elements in an array. */
  mint pOut = 0;
  for (mint i=0; i<numElements; i++) {pOut += arrIn[i];}
  return pOut;
}

mint closest_power_of_two(const mint aIn) {
  // Find the closest power of 2.
  mint lc = 0;
  while ((pow(2,lc) < aIn) && (lc < 1000000)) {lc += 1;}
  mint prev_power_of_two = pow(2, lc-1);
  mint next_power_of_two = pow(2, lc);

  mint prevDist = abs(prev_power_of_two - aIn);
  mint nextDist = abs(next_power_of_two - aIn);
  if (prevDist < nextDist) {
    return mugy_max(1,prev_power_of_two);
  } else {
    return mugy_max(1,next_power_of_two);
  }
}

mint fileExists(const char *fileNameIn) {
  /* Check if a file exists by opening it an closing it. */
  FILE *file_p;
  if ((file_p = fopen(fileNameIn, "r"))) {
    fclose(file_p);
    return true;
  }
  return false;
}

void abortSimulation(const char *errorString) {
  // Abruptly terminate the simulation (e.g. due to an error).
  printf("%s\n\n",errorString);
  exit(0);
}

// Turn elements in an array into a string, separated by a space.
void arr2str_mint(char *str, const mint *arr, const mint numElements, const char *preStr, const char *postStr) {
  str += sprintf(str, preStr);
  for (mint i=0; i<numElements; i++) str += sprintf(str, " %"fmt_mint, arr[i]);
  str += sprintf(str, postStr);
}
void arr2str_real(char *str, const real *arr, const mint numElements, const char *preStr, const char *postStr) {
  str += sprintf(str, preStr);
  for (mint i=0; i<numElements; i++) str += sprintf(str, " %"fmt_real, arr[i]);
  str += sprintf(str, postStr);
}

// Function to obtain a pointer to the i-th array in an flat array of
// multiple arrays, each array with numElem elements.
mint* getArray_mint(mint *arr, const mint *numElem, const mint i) {
  mint *newPtr = arr;
  mint off = 0;
  for (mint d=0; d<i; d++) off += numElem[d];
  return newPtr+off;
}
real* getArray_real(real *arr, const mint *numElem, const mint i) {
  real *newPtr = arr;
  mint off = 0;
  for (mint d=0; d<i; d++) off += numElem[d];
  return newPtr+off;
}
