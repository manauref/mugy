/* A mugy header file:

   Utility functions used in mugy.
*/

#ifndef UTILITIES
#define UTILITIES

#include <math.h>     // e.g. for pow.
#include <stdlib.h>   // e.g. for abs.
#include <stdbool.h>  // e.g. for bool, true, false.
#include <stdio.h>    // e.g. for fopen.
#include "data_mugy.h"
#include "macros.h"

mint prod_mint(const mint *arrIn, mint numElements);  // Product of the elements in array.
mint sum_mint(const mint *arrIn, mint numElements);  // Sum of the elements in array.

// Obtain the power of 2 closest to aIn.
mint closest_power_of_two(const mint aIn);

// Check if a file exists.
mint fileExists(const char *fileNameIn);

// Print 'errorString' and exit/terminate the simulation. 
void abortSimulation(const char *errorString);

// Turn elements in an array into a string, separated by a space.
void arr2str_mint(char *str, const mint *arr, const mint numElements, const char *preStr, const char *postStr);
void arr2str_real(char *str, const real *arr, const mint numElements, const char *preStr, const char *postStr);

// Function to obtain a pointer to the i-th array in an flat array of
// multiple arrays, each array with numElem elements.
mint* getArray_mint(mint *arr, const mint *numElem, const mint i);
real* getArray_real(real *arr, const mint *numElem, const mint i);

#endif
