/* mugy
   
   A collection of utility functions.
*/
#include "utilities.h"

int prod_int(const int *arrIn, const int numElements) {
  /* Compute the product of the elements in an array. */
  int pOut = 1;
  for (int i=0; i<numElements; i++) {pOut *= arrIn[i];}
  return pOut;
}

int closest_power_of_two(const int aIn) {
  /* Find the closest power of 2. */
  int lc = 0;
  while ((pow(2,lc) < aIn) && (lc < 1000000)) {lc += 1;}
  int prev_power_of_two = pow(2, lc-1);
  int next_power_of_two = pow(2, lc);

  int prevDist = abs(prev_power_of_two - aIn);
  int nextDist = abs(next_power_of_two - aIn);
  if (prevDist < nextDist) {
    return prev_power_of_two;
  } else {
    return next_power_of_two;
  }
}

int fileExists(const char *fileNameIn) {
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
void arr2str_int(char *str, const int *arr, const int numElements, const char *preStr, const char *postStr) {
  str += sprintf(str, preStr);
  for (int i=0; i<numElements; i++) str += sprintf(str, " %d", arr[i]);
  str += sprintf(str, postStr);
}
void arr2str_real(char *str, const real *arr, const int numElements, const char *preStr, const char *postStr) {
  str += sprintf(str, preStr);
  for (int i=0; i<numElements; i++) str += sprintf(str, " %"SCNfREAL, arr[i]);
  str += sprintf(str, postStr);
}


