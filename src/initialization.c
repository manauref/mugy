/* mugy
   
   Functions used to initialize the simulation.
*/
#include "initialization.h"

void read_inputs(const char *fileNameIn, int *NkxG, real *kxMin) {
  // Read input values from input file.
  int nFrames, mpiProcs[3];
  real kxMaxDyn[3], dt, endTime;
  FILE *file_p;

  printf("\n Reading inputs from %s\n",fileNameIn);

  file_p = fopen(fileNameIn, "r");  // Open for read only.

  fscanf(file_p, "%*s");  // spaceIn
  fscanf(file_p, "%*s %*s %d %d %d", &NkxG[0], &NkxG[1], &NkxG[2]);
#if USE_SINGLE_PRECISION > 0
  fscanf(file_p, "%*s %*s %f %f %f", &kxMin[0], &kxMin[1], &kxMin[2]);
  fscanf(file_p, "%*s %*s %f %f %f", &kxMaxDyn[0], &kxMaxDyn[1], &kxMaxDyn[2]);
  fscanf(file_p, "%*s %*s %f", &dt);
  fscanf(file_p, "%*s %*s %f", &endTime);
#else
  fscanf(file_p, "%*s %*s %lf %lf %lf", &kxMin[0], &kxMin[1], &kxMin[2]);
  fscanf(file_p, "%*s %*s %lf %lf %lf", &kxMaxDyn[0], &kxMaxDyn[1], &kxMaxDyn[2]);
  fscanf(file_p, "%*s %*s %lf", &dt);
  fscanf(file_p, "%*s %*s %lf", &endTime);
#endif
  fscanf(file_p, "%*s %*s %d", &nFrames);
  fscanf(file_p, "%*s %*s %d %d %d", &mpiProcs[0], &mpiProcs[1], &mpiProcs[2]);

  fclose(file_p);
}

void init_gridsG(const int *userNkxG, struct gridType *grid) {
  // Set number of cells in de-aliased, aliased and real space global grids.

  /* Given user-input number of distinct dealised wavenumbers, Nkx,
     the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
     to be a power of 2, so we may need to adjust Nkx. */
  printf("\n User requested  NkxG=(%d,%d,%d) distinct wavenumbers (absolute magnitude)\n", userNkxG[0], userNkxG[1], userNkxG[2]);
  grid->NxaG[0] = closest_power_of_two(3*(userNkxG[0]-1));
  grid->NxaG[1] = closest_power_of_two(3*(userNkxG[1]-1));
  grid->NxaG[2] = 1;

  // Number of distinct aliased (absolute) wavenumbers.
  grid->NkxaG[0] = grid->NxaG[0]/2+1;
  grid->NkxaG[1] = grid->NxaG[1]/2+1;
  grid->NkxaG[2] = grid->NxaG[2]/2+1;
  // Length of aliased arrays along kx and ky.
  grid->NekxaG[0] = grid->NxaG[0];
  grid->NekxaG[1] = grid->NkxaG[1];
  grid->NekxaG[2] = grid->NkxaG[2];

  // Recompute the number of distinct de-aliased (absolute) wavenumbers.
  grid->NkxG[0] = 2*(grid->NkxaG[0]-1)/3+1;
  grid->NkxG[1] = 2*(grid->NkxaG[1]-1)/3+1;
  grid->NkxG[2] = 1;
  // Length of de-aliased arrays along kx and ky.
  grid->NekxG[0] = 2*(grid->NkxG[0]-1)+1;
  grid->NekxG[1] = grid->NkxG[1];
  grid->NekxG[2] = grid->NkxG[2];

  // Number of cells in de-aliased real-space.
  grid->NxG[0] = 2*(grid->NkxG[0]-1)+1;
  grid->NxG[1] = 2*(grid->NkxG[1]-1);
  grid->NxG[2] = 1;

  printf("\n Proceeding with :\n");
  printf(" Number of distinct de-aliased absolute wavenumbers: NkxG   =%6d | NkyG   =%6d | NkzG   =%6d\n",   grid->NkxG[0],   grid->NkxG[1],   grid->NkxG[2]);
  printf(" Length of de-aliased k-space arrays:                NekxG  =%6d | NekyG  =%6d | NekzG  =%6d\n",  grid->NekxG[0],  grid->NekxG[1],  grid->NekxG[2]);
  printf(" Number of distinct aliased absolute wavenumbers:    NkxaG  =%6d | NkyaG  =%6d | NkzaG  =%6d\n",  grid->NkxaG[0],  grid->NkxaG[1],  grid->NkxaG[2]);
  printf(" Length of aliased k-space arrays:                   NekxaG =%6d | NekyaG =%6d | NekzaG =%6d\n", grid->NekxaG[0], grid->NekxaG[1], grid->NekxaG[2]);
  printf(" Number of aliased real space cells:                 NxaG   =%6d | NyaG   =%6d | NzaG   =%6d\n",   grid->NxaG[0],   grid->NxaG[1],   grid->NxaG[2]);
  printf(" Number of de-aliased real space cells:              NxG    =%6d | NyG    =%6d | NzG    =%6d\n",    grid->NxG[0],    grid->NxG[1],    grid->NxG[2]);
}
