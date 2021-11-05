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

// Product of the elements in array.
int prod_int(const int *arrIn, int numElements);

// Obtain the power of 2 closest to aIn.
int closest_power_of_two(const int aIn);

// Check if a file exists.
int fileExists(const char *fileNameIn);

// Print 'errorString' and exit/terminate the simulation. 
void abortSimulation(const char *errorString);

// Turn elements in an array into a string, separated by a space.
void arr2str_int(char *str, const int *arr, const int numElements, const char *preStr, const char *postStr);
void arr2str_real(char *str, const real *arr, const int numElements, const char *preStr, const char *postStr);

#endif
