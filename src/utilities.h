/* A mugy header file:

   Utility functions used in mugy.
*/

#ifndef UTILITIES
#define UTILITIES

#include <math.h>     // e.g. for pow.
#include <stdlib.h>   // e.g. for abs.
#include <stdbool.h>  // e.g. for bool, true, false.
#include <stdio.h>    // e.g. for fopen.

// Product of the elements in array.
int prod(const int *arrIn);

// Obtain the power of 2 closest to aIn.
int closest_power_of_two(const int aIn);

// Check if a file exists.
int fileExists(const char *fileNameIn);

// Print 'errorString' and exit/terminate the simulation. 
void abortSimulation(const char *errorString);

#endif
