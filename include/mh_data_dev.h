/* mugy: mh_data_dev.h
 *
 * Data types and some macros used in mugy device (GPU) operations.
 *
 */
#pragma once

#include "mh_macros.h"
#include <stddef.h>

// Copy mint/real/fourier data moments between host and device.
void memcpy_mint_dev(mint *dest, mint *src, mint dofs, enum memcpy_dir_dev dir);
void memcpy_real_dev(real *dest, real *src, mint dofs, enum memcpy_dir_dev dir);
void memcpy_fourier_dev(void *dest, void *src, mint dofs, enum memcpy_dir_dev dir);
void *mugy_memcpy_dev(void *dest, void *src, size_t sz, enum memcpy_dir_dev dir);
