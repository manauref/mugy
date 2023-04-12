/* mugy: mh_ffts_dev.h
 * 
 * Compute FFTs on the device.
 *
 */
#pragma once

#include "mh_ffts.h"
#include "mh_grid.h"
#include "mh_comms.h"
#include "mh_macros.h"

typedef struct mugy_fft_fam_dev mugy_fft_fam_dev;

struct mugy_fft_fam_dev *mugy_fft_init_dev(struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_population popL, struct mugy_comms comms);

void mugy_fft_c2r_dev(struct mugy_fft_fam_dev *ffts, struct mugy_array *fOut, struct mugy_array *fkIn, enum mugy_fft_type ttype);
void mugy_fft_r2c_dev(struct mugy_fft_fam_dev *ffts, struct mugy_array *fkOut, struct mugy_array *fIn, enum mugy_fft_type ttype);

void mugy_fft_terminate_dev(struct mugy_fft_fam_dev *ffts);
