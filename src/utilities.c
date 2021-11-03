/* mugy
   
   A collection of utility functions.
*/
#include "utilities.h"

int prod(const int *arrIn) {
  /* Compute the product of the elements in an array. */
  int pOut = 1;
  int num = sizeof(&arrIn)/sizeof(int);
  for (int i=0; i<num; i++) {pOut *= arrIn[i];}
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
