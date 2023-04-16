/* mugy: initialization.c
 *
 * Functions used to initialize the simulation.
 *
 */

#include <string.h>   // e.g. for strcat, strlen.
#include "mh_utilities.h"
#include "mh_alloc.h"
#include "mh_comms.h"
#include "mh_io_tools.h"
#include "mh_initialization.h"
#include "mh_initialization_dev.h"
#include "mh_data.h"
#include "mh_fourier_ho.h"
#include "mh_grid.h"

// Read an variable from input file.
void readFileVar_mint(FILE *fp, const mint numElements, mint *var) {
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%d", var);
  } else {
    for (mint i=0; i<numElements; i++) fscanf(fp, "%d", &var[i]);
  }
};
void readFileVar_real(FILE *fp, const mint numElements, real *var) {
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%"fmt_real, var);
  } else {
    for (mint i=0; i<numElements; i++) fscanf(fp, "%"fmt_real, &var[i]);
  }
};

// Read species parameter composed of numElements[s] for the s-th species.
// Skip the species prior to the s-th species, and when numElements[i] is
// greater than one allocate an array for such parameter.
void readFileSpeciesPar_mint(mint **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements) {
  fscanf(fp, "%*s = ");
  // Skip species prior to sIdx (presumably already read).
  for (mint s=0; s<sIdx; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, "%*d");
  }
  // Read this species' parameter.
  if (numElements[sIdx] == 1) {
    fscanf(fp, "%d", &(*var)[0]);
  } else {
    *var = alloc_mintArray_ho(numElements[sIdx]);
    for (mint i=0; i<numElements[sIdx]; i++) fscanf(fp, "%d", &(*var)[i]);
  }
  // Skip species after sIdx (presumably will be read later).
  for (mint s=sIdx+1; s<numSpecies; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, "%*d");
  }
}
void readFileSpeciesPar_real(real **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements) {
  fscanf(fp, "%*s = ");
  // Skip species prior to sIdx (presumably already read).
  char star_real_fmt[] = "%*";  strcat(star_real_fmt, fmt_real);
  for (mint s=0; s<sIdx; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, star_real_fmt);
  }
  // Read this species' parameter.
  if (numElements[sIdx] == 1) {
    fscanf(fp, "%"fmt_real, &(*var)[0]);
  } else {
    *var = mugy_alloc_real_ho(numElements[sIdx]);
    for (mint i=0; i<numElements[sIdx]; i++) fscanf(fp, "%"fmt_real, &(*var)[i]);
  }
  // Skip species after sIdx (presumably will be read later).
  for (mint s=sIdx+1; s<numSpecies; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, star_real_fmt);
  }
}

void read_inputFile(const char *fileNameIn, struct mugy_grid *grid, struct mugy_timeSetup *time,
                    struct mugy_population *pop, struct mugy_field *field, mint rank) {
  // Read input values from input file.

  struct mugy_pop *popG = &pop->global;

  if (rank == ioRank) {  // Only ioRank reads from input file.
    printf(" Reading inputs from %s\n\n",fileNameIn);

    FILE *file_p = fopen(fileNameIn, "r");  // Open for read only.

    fscanf(file_p, "%*s");  // &space.
    readFileVar_mint(file_p, nDim, &grid->global.deal.Nkx[0]);
    readFileVar_real(file_p, nDim, &grid->global.deal.kxMin[0]);
    readFileVar_real(file_p, nDim, grid->global.deal.kxMaxDyn);
    readFileVar_mint(file_p, nDim, grid->mpiProcs);
    fscanf(file_p, "%*s");  // /.

    fscanf(file_p, "%*s");  // &time.
    readFileVar_real(file_p, 1, &time->dt);
    readFileVar_real(file_p, 1, &time->endTime);
    readFileVar_mint(file_p, 1, &time->nFrames);
    // ARKODE parameters:
    readFileVar_mint(file_p, 1, &time->ark_kySplit);
    readFileVar_mint(file_p, 1, &time->ark_fastTableExp);
    readFileVar_mint(file_p, 1, &time->ark_fastTableImp);
    readFileVar_mint(file_p, 1, &time->ark_slowTable);
    readFileVar_real(file_p, 1, &time->ark_dtFast);
    readFileVar_real(file_p, 1, &time->ark_rtol);
    readFileVar_real(file_p, 1, &time->ark_atol);
    readFileVar_mint(file_p, 1, &time->ark_ewtScaling);
    fscanf(file_p, "%*s");  // /.

    fscanf(file_p, "%*s");  // &species.
    readFileVar_mint(file_p, 1, &popG->numSpecies);
    readFileVar_mint(file_p, 1, &pop->mpiProcs);
    popG->pars = (struct mugy_species_pars*) calloc(popG->numSpecies, sizeof(struct mugy_species_pars));
    mint *specNumMoms = (mint*) calloc(popG->numSpecies, sizeof(mint));
    mint *specOnes    = (mint*) calloc(popG->numSpecies, sizeof(mint));
    mint *specnDim    = (mint*) calloc(popG->numSpecies, sizeof(mint));
    readFileVar_mint(file_p, popG->numSpecies, specNumMoms);
    mint filePos = ftell(file_p);
    for (mint s=0; s<popG->numSpecies; s++) {
      popG->pars[s].numMoments = specNumMoms[s];
      specOnes[s] = 1;
      specnDim[s] = nDim;
    }
    real* real_p; real** real_pp; mint* mint_p;
    for (mint s=0; s<popG->numSpecies; s++) {
      fseek(file_p, filePos, SEEK_SET);
      real_p  =    &popG->pars[s].qCharge; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =     &popG->pars[s].muMass; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =        &popG->pars[s].tau; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =       &popG->pars[s].omSt; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =        &popG->pars[s].omd; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =      &popG->pars[s].delta; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =  &popG->pars[s].deltaPerp; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =        &popG->pars[s].eta; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_pp =      &popG->pars[s].alpha; readFileSpeciesPar_real(real_pp, file_p, s, popG->numSpecies, specNumMoms);
      real_pp =         &popG->pars[s].nu; readFileSpeciesPar_real(real_pp, file_p, s, popG->numSpecies, specNumMoms);
      real_p  =     &popG->pars[s].delta0; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_pp = &popG->pars[s].hDiffOrder; readFileSpeciesPar_real(real_pp, file_p, s, popG->numSpecies, specnDim   );
      real_pp =      &popG->pars[s].hDiff; readFileSpeciesPar_real(real_pp, file_p, s, popG->numSpecies, specnDim   );
      real_pp =   &popG->pars[s].kDiffMin; readFileSpeciesPar_real(real_pp, file_p, s, popG->numSpecies, specnDim   );
      mint_p  =       &popG->pars[s].icOp; readFileSpeciesPar_mint(&mint_p, file_p, s, popG->numSpecies, specOnes   );
      real_pp =    &popG->pars[s].initAux; readFileSpeciesPar_real(real_pp, file_p, s, popG->numSpecies, specnDim   );
      real_p  =      &popG->pars[s].initA; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
      real_p  =     &popG->pars[s].noiseA; readFileSpeciesPar_real(&real_p, file_p, s, popG->numSpecies, specOnes   );
    }
    free(specnDim); free(specOnes); free(specNumMoms);
    fscanf(file_p, "%*s");  // /.
  
    fscanf(file_p, "%*s");  // &field.
    readFileVar_real(file_p, 1, &field->pars.lambdaD);
    readFileVar_mint(file_p, 1, &field->pars.pade);
    readFileVar_mint(file_p, 1, &field->pars.icOp);
    fscanf(file_p, "%*s");  // /.

    fclose(file_p);
  }

  // Broadcast to other processes.
  MPI_Bcast(&grid->global.deal.Nkx     , nDim, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->global.deal.kxMin   , nDim, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->global.deal.kxMaxDyn, nDim, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->mpiProcs   , nDim, mpi_mint, ioRank, MPI_COMM_WORLD);

  MPI_Bcast(&time->dt              , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->endTime         , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->nFrames         , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_kySplit     , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_fastTableExp, 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_fastTableImp, 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_slowTable   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_dtFast      , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_rtol        , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_atol        , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_ewtScaling  , 1, mpi_mint, ioRank, MPI_COMM_WORLD);

  MPI_Bcast(&popG->numSpecies, 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&pop->mpiProcs   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  if (rank != ioRank) popG->pars = (struct mugy_species_pars*) calloc(popG->numSpecies, sizeof(struct mugy_species_pars));
  for (mint s=0; s<popG->numSpecies; s++) {
    MPI_Bcast(&popG->pars[s].numMoments,                       1, mpi_mint, ioRank, MPI_COMM_WORLD);
    if (rank != ioRank) {
      popG->pars[s].alpha      = mugy_alloc_real_ho(popG->pars[s].numMoments);
      popG->pars[s].nu         = mugy_alloc_real_ho(popG->pars[s].numMoments);
      popG->pars[s].hDiffOrder = mugy_alloc_real_ho(nDim);
      popG->pars[s].hDiff      = mugy_alloc_real_ho(nDim);
      popG->pars[s].kDiffMin   = mugy_alloc_real_ho(nDim);
      popG->pars[s].initAux    = mugy_alloc_real_ho(nDim);
    }
    MPI_Bcast(&popG->pars[s].qCharge   ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].muMass    ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].tau       ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].omSt      ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].omd       ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].delta     ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].deltaPerp ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].eta       ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(popG->pars[s].alpha      , popG->pars[s].numMoments, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(popG->pars[s].nu         , popG->pars[s].numMoments, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].delta0    ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(popG->pars[s].hDiffOrder ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(popG->pars[s].hDiff      ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(popG->pars[s].kDiffMin   ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].icOp      ,                       1, mpi_mint, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(popG->pars[s].initAux    ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].initA     ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&popG->pars[s].noiseA    ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
  }

  MPI_Bcast(&field->pars.lambdaD, 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&field->pars.pade   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&field->pars.icOp   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);

}

void read_inputs(mint argc, char *argv[], struct mugy_ioSetup *ioSet, struct mugy_grid *grid, struct mugy_timeSetup *time,
                 struct mugy_population *pop, struct mugy_field *field, mint rank) {
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

  read_inputFile(ioSet->inputFile, grid, time, pop, field, rank);

  // Set the total number of moments.
  struct mugy_pop *popG = &pop->global;
  popG->numMomentsTot = 0;
  for (mint s=0; s<popG->numSpecies; s++) popG->numMomentsTot += popG->pars[s].numMoments;

}

// Initialize device
void device_init(struct mugy_comms *comms) {
  device_init_dev(comms);
}

void set_initialConditions(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid,
  struct mugy_ffts *fftMan, struct mugy_ioManager *ioman) {
  // Impose the initial conditions on the moments and thoe potential.

  struct mugy_array *momk = pop->local.momk[0]; // Put ICs in first stepper field.

  // NOTE: For now assume initialOp is the same for all species.
  mint initialOp = pop->local.pars[0].icOp; 

  if (initialOp == 0) {
    // Initialize in real space and transform to Fourier.
    struct mugy_realGrid *gridL = &grid->local.deal.dual;
    struct mugy_array *momIC = mugy_population_alloc_realMoments(*gridL, pop->local, hostAndDeviceMem);

    for (mint s=0; s<pop->local.numSpecies; s++) {
      real initA = pop->local.pars[s].initA;

      real *den_p  = mugy_population_getMoment_real(*gridL, pop->local, s, denIdx, momIC->ho);  // Get density of species s.
      real *temp_p = mugy_population_getMoment_real(*gridL, pop->local, s, tempIdx, momIC->ho);  // Get temperature of species s.

      for (mint linIdx=0; linIdx<gridL->NxTot; linIdx++) {
        mint xIdx[nDim];
        mugy_grid_lin2sub_real(&xIdx[0], linIdx, *gridL);  // Convert linear index to multidimensional x index.
        real x[nDim];
        mugy_grid_get_x(&x[0], xIdx, *gridL);

        // Initial density: a superposition of sines and cosines.
        double kx = grid->local.deal.kxMin[0];
        double ky = grid->local.deal.kxMin[1];
        den_p[0] += initA*sin(kx*x[0])*cos(ky*x[1]);
        den_p++;
        
        // Initial temperature = 0.
        temp_p[0] = 0.;
        temp_p++;
      }
    }

    // Copy initialized moments from host to device.
    mugy_array_copy(momIC, momIC, host2device);

    // Forward FFT moments.
    mugy_fft_r2c(fftMan, momk, momIC, mugy_fft_mom_xy, deviceComp);

//    //......................................................
//    // Test FFT of a moments.
//    struct mugy_array *momICk = mugy_population_alloc_fourierMoments(grid->local.deal, pop->local, hostAndDeviceMem);
//
////    mugy_fft_r2c(fftMan, momICk, momIC, mugy_fft_mom_xy, hostComp);
////    mugy_fft_c2r(fftMan, momIC, momICk, mugy_fft_mom_xy, hostComp);
//    mugy_fft_r2c(fftMan, momICk, momIC, mugy_fft_mom_xy, deviceComp);
//    mugy_fft_c2r(fftMan, momIC, momICk, mugy_fft_mom_xy, deviceComp);
//    mugy_array_copy(momIC, momIC, device2host);
//
//    struct mugy_ad_file *fh = mugy_io_create_moments_file(ioman, "mom", *grid, *pop, real_enum);
//    mugy_io_write_mugy_array(NULL, "mom", fh, momIC);
//    mugy_io_close_file(fh);
//
//    mugy_array_free(momICk, hostAndDeviceMem);
//    //......................................................
//
//    //......................................................
//    // Test FFT of a single array
//    struct mugy_array *fxy_r = mugy_array_alloc(real_enum, gridL->NxTot, hostAndDeviceMem);
//    struct mugy_array *fxy_k = mugy_array_alloc(fourier_enum, grid->local.deal.NekxTot, hostAndDeviceMem);
//
//    // Assign real array to a linear combo of sines and cosines.
//    real *fxy_rp = fxy_r->ho;
//    for (mint linIdx=0; linIdx<gridL->NxTot; linIdx++) {
//      real initA = pop->local.pars[0].initA;
//      mint xIdx[nDim];
//      mugy_grid_lin2sub_real(&xIdx[0], linIdx, *gridL);  // Convert linear index to multidimensional x index.
//      real x[nDim];
//      mugy_grid_get_x(&x[0], xIdx, *gridL);
//
//      fxy_rp[0] = 0.;
//      double kx = grid->local.deal.kxMin[0];
//      double ky = grid->local.deal.kxMin[1];
//      fxy_rp[0] += initA*sin(kx*x[0])*cos(ky*x[1]);
//      fxy_rp++;
//    }
//
////    mugy_fft_r2c(fftMan, fxy_k, fxy_r, mugy_fft_xy, hostComp);
////    mugy_fft_c2r(fftMan, fxy_r, fxy_k, mugy_fft_xy, hostComp);
//
//    mugy_array_copy(fxy_r, fxy_r, host2device);
//    mugy_fft_r2c(fftMan, fxy_k, fxy_r, mugy_fft_xy, deviceComp);
//    mugy_fft_c2r(fftMan, fxy_r, fxy_k, mugy_fft_xy, deviceComp);
//    mugy_array_copy(fxy_r, fxy_r, device2host);
//
//    struct mugy_ad_file *fhr = mugy_io_create_mugy_array_file(ioman, "arr", gridL, real_enum);
//    mugy_io_write_mugy_array(NULL, "arr", fhr, fxy_r);
//    mugy_io_close_file(fhr);
//
//    mugy_array_free(fxy_r, hostAndDeviceMem);
//    mugy_array_free(fxy_k, hostAndDeviceMem);
//    //......................................................

    mugy_array_free(momIC, hostAndDeviceMem);

  } else if (initialOp == 1) {
    // Initialize with a k-spce power law.
    struct mugy_fourierGrid *gridL = &grid->local.deal;
    real *kxMin = &gridL->kxMin[0];

    for (mint s=0; s<pop->local.numSpecies; s++) {
      real initA    = pop->local.pars[s].initA;
      real *initAux = &pop->local.pars[s].initAux[0];

      // Get density and temperature of species s.
      fourier *den_p  = mugy_population_getMoment_fourier(*gridL, pop->local, s, denIdx, momk->ho);
      fourier *temp_p = mugy_population_getMoment_fourier(*gridL, pop->local, s, tempIdx, momk->ho);

      for (mint linIdx=0; linIdx<gridL->NekxTot; linIdx++) {
        mint kxIdx[nDim];
        mugy_grid_lin2sub_fourier(&kxIdx[0], linIdx, *gridL);  // Convert linear index to multidimensional kx index.
        real kx[nDim];
        mugy_grid_get_kx(&kx[0], kxIdx, *gridL);
  
        // Set density to a power-law in k-space.
        den_p[0] = initA*(pow((kxMin[0]+fabs(kx[0]))/kxMin[0],initAux[0]))
                        *(pow((kxMin[1]+fabs(kx[1]))/kxMin[1],initAux[1]));
        den_p++;
  
        // Set the initial temperature (fluctuations) to zero.
        temp_p[0] = 0.;
        temp_p++;
      };
    };

    // Copy initialized moments from host to device.
    mugy_array_copy(momk, momk, host2device);
  }

  // Solve the Poisson equation to compute phik.
//  mugy_array_copy(momk, momk, device2host);
  mugy_field_poisson_solve(field, pop, grid, 0);

  mugy_array_copy(field->phik, field->phik, device2host);
  mugy_io_write_mugy_array(ioman, "phik", NULL, field->phik);

}
