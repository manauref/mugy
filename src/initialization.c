/* mugy
   
   Functions used to initialize the simulation.
*/
#include "initialization.h"

void readFileVar_real(FILE *fp, const int numElements, real *var) {
  // Read a real variable from input file.
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%"SCNfREAL, var);
  } else {
    for (int i=0; i<numElements; i++) fscanf(fp, "%"SCNfREAL, &var[i]);
  }
};
void readFileVar_int(FILE *fp, const int numElements, int *var) {
  // Read an int variable from input file.
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%d", var);
  } else {
    for (int i=0; i<numElements; i++) fscanf(fp, "%d", &var[i]);
  }
};
void allocAndReadFileVar_real(FILE *fp, const int numElements, real **var) {
  // Allocate a real array and read its elements from input file.
  fscanf(fp, "%*s = ");
  *var = alloc_realArray(numElements);
  for (int i=0; i<numElements; i++) fscanf(fp, "%"SCNfREAL, &(*var)[i]);
};
void allocAndReadFileVar_int(FILE *fp, const int numElements, int **var) {
  // Allocate an int array and read its elements from input file.
  fscanf(fp, "%*s = ");
  *var = alloc_intArray(numElements);
  for (int i=0; i<numElements; i++) fscanf(fp, "%d", &(*var)[i]);
};

void read_inputFile(const char *fileNameIn, struct grid *grid, struct timeSetup *time,
                    struct speciesParameters *spec, struct fieldParameters *field) {
  // Read input values from input file.

  if (myRank == 0) printf(" Reading inputs from %s\n\n",fileNameIn);

  FILE *file_p = fopen(fileNameIn, "r");  // Open for read only.

  fscanf(file_p, "%*s");  // &space.
  readFileVar_int( file_p, 3, grid->fG.Nkx);
  readFileVar_real(file_p, 3, grid->fG.kxMin);
  readFileVar_real(file_p, 3, grid->fG.kxMaxDyn);
  readFileVar_int( file_p, 3, grid->mpiProcs);
  fscanf(file_p, "%*s");  // /.

  fscanf(file_p, "%*s");  // &time.
  readFileVar_real(file_p, 1, &time->dt);
  readFileVar_real(file_p, 1, &time->endTime);
  readFileVar_int( file_p, 1, &time->nFrames);
  // ARKODE parameters:
  readFileVar_int( file_p, 1, &time->ark_kySplit);
  readFileVar_int( file_p, 1, &time->ark_fastTableExp);
  readFileVar_int( file_p, 1, &time->ark_fastTableImp);
  readFileVar_int( file_p, 1, &time->ark_slowTable);
  readFileVar_real(file_p, 1, &time->ark_dtFast);
  readFileVar_real(file_p, 1, &time->ark_rtol);
  readFileVar_real(file_p, 1, &time->ark_atol);
  readFileVar_int( file_p, 1, &time->ark_ewtScaling);
  fscanf(file_p, "%*s");  // /.

  fscanf(file_p, "%*s");  // &species.
  readFileVar_int( file_p, 1, &spec->numSpecies);
  readFileVar_int( file_p, 1, &spec->numMoments);
  readFileVar_int( file_p, 1, &spec->mpiProcs);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,    &spec->qCharge);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,     &spec->muMass);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,        &spec->tau);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,       &spec->omSt);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,        &spec->omd);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,      &spec->delta);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,  &spec->deltaPerp);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,        &spec->eta);
  allocAndReadFileVar_real(file_p, spec->numMoments*spec->numSpecies,      &spec->alpha);
  allocAndReadFileVar_real(file_p, spec->numMoments*spec->numSpecies,         &spec->nu);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,     &spec->delta0);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies, &spec->hDiffOrder);
  allocAndReadFileVar_real(file_p,                3*spec->numSpecies,      &spec->hDiff);
  allocAndReadFileVar_real(file_p,                3*spec->numSpecies,   &spec->kDiffMin);
  allocAndReadFileVar_int( file_p,                  spec->numSpecies,       &spec->icOp);
  allocAndReadFileVar_real(file_p,                3*spec->numSpecies,    &spec->initAux);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,      &spec->initA);
  allocAndReadFileVar_real(file_p,                  spec->numSpecies,     &spec->noiseA);
  fscanf(file_p, "%*s");  // /.

  fscanf(file_p, "%*s");  // &field.
  readFileVar_real(file_p, 1, &field->lambdaD);
  readFileVar_int( file_p, 1, &field->pade);
  readFileVar_int( file_p, 1, &field->icOp);
  fscanf(file_p, "%*s");  // /.

  fclose(file_p);

}

void read_inputs(int argc, char *argv[], struct ioSetup *ioSet, struct grid *grid, struct timeSetup *time,
                 struct speciesParameters *spec, struct fieldParameters *field) {
  // Read inputs from command line arguments and input file.

  // Check for commandline arguments.
  // Currently we expect two arguments in this order:
  //   1) name of the input file.
  //   2) absolute address of the output directory.
  // For restarts add the restart directory as a 3rd argument.
  if (argc < 3) {
    printf("\n --> Not enough inputs. Need input file and output directory as command line arguments.\n");
    abortSimulation(" TERMINATING ");
  } else {
    ioSet->inputFile = argv[1];
    ioSet->outputDir = argv[2];
    ioSet->isRestart   = false;
    ioSet->outToOldDir = false;
    if (argc > 3) {  // Simulation is a restart of a previous one.
      ioSet->restartDir = argv[3];
      ioSet->isRestart  = true;
      // Output to the same directory as the previous run?
      char *checkFile = malloc(strlen(ioSet->outputDir)+strlen("phik.bp"));
      checkFile[0] = '\0';
      strcat(strcat(checkFile,ioSet->outputDir),"phik.bp");
      ioSet->outToOldDir = fileExists(checkFile);
      free(checkFile);
    }
  }

  read_inputFile(ioSet->inputFile, grid, time, spec, field);

}

void init_grids(struct grid *grid) {
  // Set number of cells in de-aliased, aliased and real space global grids.

  /* Given user-input number of distinct dealised wavenumbers, Nkx,
     the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
     to be a power of 2, so we may need to adjust Nkx. */
  arrPrint_int(grid->fG.Nkx, nDim, " User requested  NkxG=(",") distinct wavenumbers (absolute magnitude)\n");
  grid->fGa.dual.Nx[0] = closest_power_of_two(3*(grid->fG.Nkx[0]-1));
  grid->fGa.dual.Nx[1] = closest_power_of_two(3*(grid->fG.Nkx[1]-1));
  grid->fGa.dual.Nx[2] = 1;

  // Number of distinct aliased (absolute) wavenumbers.
  grid->fGa.Nkx[0] = grid->fGa.dual.Nx[0]/2+1;
  grid->fGa.Nkx[1] = grid->fGa.dual.Nx[1]/2+1;
  grid->fGa.Nkx[2] = grid->fGa.dual.Nx[2]/2+1;
  // Length of aliased arrays along kx and ky.
  grid->fGa.Nekx[0] = grid->fGa.dual.Nx[0];
  grid->fGa.Nekx[1] = grid->fGa.Nkx[1];
  grid->fGa.Nekx[2] = grid->fGa.Nkx[2];

  // Recompute the number of distinct de-aliased (absolute) wavenumbers.
  grid->fG.Nkx[0] = 2*(grid->fGa.Nkx[0]-1)/3+1;
  grid->fG.Nkx[1] = 2*(grid->fGa.Nkx[1]-1)/3+1;
  grid->fG.Nkx[2] = 1;
  // Length of de-aliased arrays along kx and ky.
  grid->fG.Nekx[0] = 2*(grid->fG.Nkx[0]-1)+1;
  grid->fG.Nekx[1] = grid->fG.Nkx[1];
  grid->fG.Nekx[2] = grid->fG.Nkx[2];

  // Number of cells in de-aliased real-space.
  grid->fG.dual.Nx[0] = 2*(grid->fG.Nkx[0]-1)+1;
  grid->fG.dual.Nx[1] = 2*(grid->fG.Nkx[1]-1);
  grid->fG.dual.Nx[2] = 1;

  r0printf("\n Proceeding with :\n");
  arrPrint_int(grid->fG.Nkx,      nDim, " Number of distinct de-aliased absolute wavenumbers: NkxG   =", "\n");
  arrPrint_int(grid->fG.Nekx,     nDim, " Length of de-aliased k-space arrays:                NekxG  =", "\n");
  arrPrint_int(grid->fGa.Nkx,     nDim, " Number of distinct aliased absolute wavenumbers:    NkxaG  =", "\n");
  arrPrint_int(grid->fGa.Nekx,    nDim, " Length of aliased k-space arrays:                   NekxaG =", "\n");
  arrPrint_int(grid->fGa.dual.Nx, nDim, " Number of aliased real space cells:                 NxaG   =", "\n");
  arrPrint_int(grid->fG.dual.Nx,  nDim, " Number of de-aliased real space cells:              NxG    =", "\n");

  arrPrint_real(grid->fG.kxMin,    nDim, " Minimum absolute magnitude of wavenumbers: kxMin    =", "\n");
  arrPrint_real(grid->fG.kxMaxDyn, nDim, " Largest wavenumbers evolved:               kxMaxDyn =", "\n");

}

void free_speciesPars(struct speciesParameters *spec) {
  // Deallocate memory used in species struct.
  free(  spec->qCharge);
  free(   spec->muMass);
  free(      spec->tau);
  free(     spec->omSt);
  free(      spec->omd);
  free(    spec->delta);
  free(spec->deltaPerp);
  free(      spec->eta);
  free(    spec->alpha);
  free(       spec->nu);
  free(   spec->delta0);
  free(spec->hDiffOrder);
  free(    spec->hDiff);
  free( spec->kDiffMin);
  free(     spec->icOp);
  free(  spec->initAux);
  free(    spec->initA);
  free(   spec->noiseA);
}
