/* mugy: mh_data.h
 *
 * Methods to manipulate generic data.
 *
 */
#pragma once

#include "mh_macros.h"
#include <stddef.h>

// Copy mint/real/fourier data (between host and device, or within a host or device).
void memcpy_mint(mint *dest, mint *src, mint numElements, enum memcpy_dir_dev dir);
void memcpy_real(real *dest, real *src, mint numElements, enum memcpy_dir_dev dir);
void memcpy_fourier(void *dest, void *src, mint numElements, enum memcpy_dir_dev dir);
void *mugy_memcpy(void *dest, void *src, size_t sz, enum memcpy_dir_dev dir);
