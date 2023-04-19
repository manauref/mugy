/* mugy: mpi_comms

   Communication-related operations.
*/

#include "mh_comms.h"
#include <mpi.h>
#include "mh_utilities.h"
#include "mh_alloc.h"
#include "mh_io_utilities.h"
#include <string.h>   // e.g. for memcpy.

struct mugy_comms *mugy_comms_init(mint argc, char *argv[]) {
  // Initialize MPI, get rank of this process and total number of processes.
  struct mugy_comms *comms = (struct mugy_comms *) malloc(sizeof(struct mugy_comms));
  comms->world = (struct mugy_comms_sub *) malloc(sizeof(struct mugy_comms_sub));
  comms->world->comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);  // Initialize MPI.
  MPI_Comm_rank(comms->world->comm, &comms->world->rank);  // Rank of this MPI process.
  MPI_Comm_size(comms->world->comm, &comms->world->size);  // Number of MPI processes.
  return comms;
}

void mugy_comms_sub_init(struct mugy_comms *comms, struct mugy_grid *grid, struct mugy_population *pop) {
  // Initialize the various sub-communicators needed.
  
  // Check the number of MPI processes is correct.
  if (prod_mint(grid->mpiProcs,nDim)*pop->mpiProcs != comms->world->size) {
    printf(" Number of MPI processes in input file (%d) differs from that in mpirun (%d).\n",
           prod_mint(grid->mpiProcs,nDim)*pop->mpiProcs, comms->world->size);
    abortSimulation(" Terminating...\n");
  }

  struct mugy_comms_sub *scomm;  // Temp pointer for ease of notation.

  // Finish initializing the world subcomm (coord initialized further down).
  scomm         = comms->world; 
  scomm->dim    = nDim+1;
  scomm->coord  = alloc_mintArray_ho(scomm->dim);
  scomm->decomp = alloc_mintArray_ho(scomm->dim);
  for (mint d=0; d<nDim; d++) scomm->decomp[d] = grid->mpiProcs[d];
  scomm->decomp[nDim] = pop->mpiProcs;
  arrPrintS_mint(scomm->decomp,nDim+1, " MPI processes along X,Y,Z,s: ", "\n", comms->world->rank);
  r0printf("\n", comms->world->rank);

  // Create a 4D communicator (3D space + species).
  mint commOrg[nDim+1] = {2,3,1,0};  // Organize the dimensions as s,Z,X,Y.
  mint bc[nDim+1]      = {true, true, true, false};  // Boundary conditions (true=periodic).
  mint reorder         = true;  // Let MPI assign arbitrary ranks.
  comms->sub4d  = (struct mugy_comms_sub *) malloc(sizeof(struct mugy_comms_sub));
  scomm         = &comms->sub4d[0]; 
  scomm->dim    = nDim+1;
  scomm->coord  = alloc_mintArray_ho(scomm->dim);
  scomm->decomp = alloc_mintArray_ho(scomm->dim);
  mint decompOrg[nDim+1], bcOrg[nDim+1];
  for (mint d=0; d<scomm->dim; d++) {
    scomm->decomp[d]      = comms->world->decomp[d];
    decompOrg[commOrg[d]] = comms->world->decomp[d];
    bcOrg[commOrg[d]]     = bc[d];
  }
  MPI_Cart_create(comms->world->comm, scomm->dim, decompOrg, bcOrg, reorder, (MPI_Comm*) &scomm->comm);
  MPI_Comm_size(scomm->comm, &scomm->size);
  MPI_Comm_rank(scomm->comm, &scomm->rank);
  MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  memcpy_mint(comms->world->coord, scomm->coord, nDim+1, MUGY_HOST2HOST);  // Make coords the same in world->

  // 3D subcommunicator (xyz space, no need for xys, xzs, yzs comms).
  comms->sub3d = (struct mugy_comms_sub *) malloc(sizeof(struct mugy_comms_sub));
  for (mint d=nDim; d<nDim+1; d++) {
    mint remain[nDim+1] = {true,true,true,true};
    remain[commOrg[d]] = false;

    scomm         = &comms->sub3d[0]; 
    scomm->dim    = 3;
    scomm->coord  = alloc_mintArray_ho(scomm->dim);
    scomm->decomp = alloc_mintArray_ho(scomm->dim);
    for (mint e=0; e<scomm->dim; e++)
      scomm->decomp[e] = comms->world->decomp[e];
    MPI_Cart_sub(comms->sub4d[0].comm, remain, (MPI_Comm*) &scomm->comm);
    MPI_Comm_size(scomm->comm, &scomm->size);
    MPI_Comm_rank(scomm->comm, &scomm->rank);
    MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  }

  // 2D subcommunicators (no need for xs,ys,zs subcomms at this point).
  comms->sub2d = (struct mugy_comms_sub *) malloc(nDim*sizeof(struct mugy_comms_sub));
  for (mint d=nDim-1; d>-1; d--) {
    mint remain[nDim+1] = {true,true,true,true};
    remain[commOrg[d]] = false;  remain[commOrg[nDim]] = false;

    scomm         = &comms->sub2d[d]; 
    scomm->dim    = 2;
    scomm->coord  = alloc_mintArray_ho(scomm->dim);
    scomm->decomp = alloc_mintArray_ho(scomm->dim);
    for (mint e=0; e<scomm->dim; e++)
      scomm->decomp[e] = comms->world->decomp[e];
    MPI_Cart_sub(comms->sub4d[0].comm, remain, (MPI_Comm*) &scomm->comm);
    MPI_Comm_size(scomm->comm, &scomm->size);
    MPI_Comm_rank(scomm->comm, &scomm->rank);
    MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  }

  // 1D subcommunicators. Organize them in s,Z,X,Y, order so processes of a
  // single X-Y plane of a single species are contiguous.
  comms->sub1d = (struct mugy_comms_sub *) malloc((nDim+1)*sizeof(struct mugy_comms_sub));
  for (mint d=0; d<nDim+1; d++) {
    mint remain[nDim+1] = {false,false,false,false};
    remain[commOrg[d]] = true;

    scomm            = &comms->sub1d[d]; 
    scomm->dim       = 1;
    scomm->coord     = alloc_mintArray_ho(scomm->dim);
    scomm->decomp    = alloc_mintArray_ho(scomm->dim);
    scomm->decomp[0] = comms->world->decomp[d];
    MPI_Cart_sub(comms->sub4d[0].comm, remain, (MPI_Comm*) &scomm->comm);
    MPI_Comm_size(scomm->comm, &scomm->size);
    MPI_Comm_rank(scomm->comm, &scomm->rank);
    MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  }

}

void distribute1dDOFs(const mint procs, const mint procID, const mint globalDOFs, mint *localDOFs, mint *firstDOF) {
  /* Distribute 1D degrees of freedom amongst processes. Give every process at
     least the floor of the division of the degrees of freedom by the number
     of processors. Then distribute the remaining elements one at a time.
       procs:  number of processes
       procID: this process' rank.
       globalDOFs: number of degrees of freedom in global domain.
       localDOFs: number of degrees of freedom in this process.
       firstDOF: first degree of freedom from the global DOFs in this process.
   */
  *localDOFs = floor((real)(globalDOFs)/(real)(procs));
  mint remainingDOFs = globalDOFs - procs*(*localDOFs);
  *firstDOF = procID*(*localDOFs);

  while (remainingDOFs > 0) {
    for (mint p=0; p<procs; p++) {
      remainingDOFs += -1;
      if (p == procID) {
        *localDOFs += 1;
      } else if (p < procID) {
        *firstDOF += 1;
      }
      if (remainingDOFs == 0) break;
    }
  };

}

void mugy_comms_distributeDOFs(struct mugy_comms *comms, struct mugy_grid *grid, struct mugy_population *pop) {
  // Distribute s,Z,X,Y amongst MPI processes.
  
  // Distribute the species.
  mint rank_s = comms->sub1d[nDim].rank;
  struct mugy_population_species *popG = pop->global,  *popL = pop->local;;
  distribute1dDOFs(pop->mpiProcs, rank_s, popG->numSpecies, &popL->numSpecies, &popL->globalSpecOff);
  popL->globalMomOff = 0;
  for (mint s=0; s<popL->globalSpecOff; s++) popL->globalMomOff += popG->pars[s].numMoments;
  popL->pars = (struct mugy_population_species_pars*) calloc(popL->numSpecies, sizeof(struct mugy_population_species_pars));
  for (mint s=0; s<popL->numSpecies; s++) {
    popL->pars[s].numMoments = popG->pars[s+popL->globalSpecOff].numMoments;

    popL->pars[s].alpha      = mugy_alloc_real_ho(popL->pars[s].numMoments);
    popL->pars[s].nu         = mugy_alloc_real_ho(popL->pars[s].numMoments);
    popL->pars[s].hDiffOrder = mugy_alloc_real_ho(nDim);
    popL->pars[s].hDiff      = mugy_alloc_real_ho(nDim);
    popL->pars[s].kDiffMin   = mugy_alloc_real_ho(nDim);
    popL->pars[s].initAux    = mugy_alloc_real_ho(nDim);

    popL->pars[s].qCharge    = popG->pars[s+popL->globalSpecOff].qCharge   ;
    popL->pars[s].muMass     = popG->pars[s+popL->globalSpecOff].muMass    ;
    popL->pars[s].tau        = popG->pars[s+popL->globalSpecOff].tau       ;
    popL->pars[s].omSt       = popG->pars[s+popL->globalSpecOff].omSt      ;
    popL->pars[s].omd        = popG->pars[s+popL->globalSpecOff].omd       ;
    popL->pars[s].delta      = popG->pars[s+popL->globalSpecOff].delta     ;
    popL->pars[s].deltaPerp  = popG->pars[s+popL->globalSpecOff].deltaPerp ;
    popL->pars[s].eta        = popG->pars[s+popL->globalSpecOff].eta       ;
    memcpy(popL->pars[s].alpha, popG->pars[s+popL->globalSpecOff].alpha, popL->pars[s].numMoments*sizeof(real));
    memcpy(popL->pars[s].nu   , popG->pars[s+popL->globalSpecOff].nu   , popL->pars[s].numMoments*sizeof(real));
    popL->pars[s].delta0     = popG->pars[s+popL->globalSpecOff].delta0    ;
    memcpy(popL->pars[s].hDiffOrder, popG->pars[s+popL->globalSpecOff].hDiffOrder, nDim*sizeof(real));
    memcpy(popL->pars[s].hDiff     , popG->pars[s+popL->globalSpecOff].hDiff     , nDim*sizeof(real));
    memcpy(popL->pars[s].kDiffMin  , popG->pars[s+popL->globalSpecOff].kDiffMin  , nDim*sizeof(real));
    popL->pars[s].icOp       = popG->pars[s+popL->globalSpecOff].icOp      ;
    memcpy(popL->pars[s].initAux, popG->pars[s+popL->globalSpecOff].initAux, nDim*sizeof(real));
    popL->pars[s].initA      = popG->pars[s+popL->globalSpecOff].initA     ;
    popL->pars[s].noiseA     = popG->pars[s+popL->globalSpecOff].noiseA    ;
  }
  // Set the total number of moments.
  popL->numMomentsTot = 0;
  for (int s=0; s<popL->numSpecies; s++) popL->numMomentsTot += popL->pars[s].numMoments;

  // Distribute the real-space and Fourier-space points. 
  struct mugy_grid_chart *gridG = grid->global,  *gridL = grid->local;
  for (mint d=0; d<nDim; d++) {
    mint rank_i = comms->sub1d[d].rank;
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->fourier->Nx[d]  , &gridL->fourier->Nx[d]  , &gridL->fourier->globalOff[d]);
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->real->Nx[d]     , &gridL->real->Nx[d]     , &gridL->real->globalOff[d]);
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->fourierAl->Nx[d], &gridL->fourierAl->Nx[d], &gridL->fourierAl->globalOff[d]);
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->realAl->Nx[d]   , &gridL->realAl->Nx[d]   , &gridL->realAl->globalOff[d]);
  }
  gridL->fourier->NxTot    = prod_mint(gridL->fourier->Nx,nDim);
  gridL->fourier->NxyTot   = prod_mint(gridL->fourier->Nx,2);
  gridL->real->NxTot       = prod_mint(gridL->real->Nx,nDim);
  gridL->real->NxyTot      = prod_mint(gridL->real->Nx,2);
  gridL->fourierAl->NxTot  = prod_mint(gridL->fourierAl->Nx,nDim);
  gridL->fourierAl->NxyTot = prod_mint(gridL->fourierAl->Nx,2);
  gridL->realAl->NxTot     = prod_mint(gridL->realAl->Nx,nDim);
  gridL->realAl->NxyTot    = prod_mint(gridL->realAl->Nx,2);

  // Create local real-space and Fourier-space coordinate arrays.
  gridL->fourier->x   = mugy_alloc_real_ho(sum_mint(gridL->fourier->Nx, nDim));
  gridL->real->x      = mugy_alloc_real_ho(sum_mint(gridL->real->Nx, nDim));
  gridL->fourierAl->x = mugy_alloc_real_ho(sum_mint(gridL->fourierAl->Nx, nDim));
  gridL->realAl->x    = mugy_alloc_real_ho(sum_mint(gridL->realAl->Nx, nDim));
  for (mint d=0; d<nDim; d++) {
    memcpy(getArray_real(gridL->fourier->x,gridL->fourier->Nx,d),
           getArray_real(gridG->fourier->x,gridG->fourier->Nx,d)+gridL->fourier->globalOff[d], gridL->fourier->Nx[d]*sizeof(real));
    memcpy(getArray_real(gridL->real->x,gridL->real->Nx,d),
           getArray_real(gridG->real->x,gridG->real->Nx,d)+gridL->real->globalOff[d], gridL->real->Nx[d]*sizeof(real));
    memcpy(getArray_real(gridL->fourierAl->x,gridL->fourierAl->Nx,d),
           getArray_real(gridG->fourierAl->x,gridG->fourierAl->Nx,d)+gridL->fourierAl->globalOff[d], gridL->fourierAl->Nx[d]*sizeof(real));
    memcpy(getArray_real(gridL->realAl->x,gridL->realAl->Nx,d),
           getArray_real(gridG->realAl->x,gridG->realAl->Nx,d)+gridL->realAl->globalOff[d], gridL->realAl->Nx[d]*sizeof(real));
  }

  // Copy global constants into local grids.
  for (mint d=0; d<nDim; d++) {
    gridL->fourier->dx[d]   = gridG->fourier->dx[d];
    gridL->real->dx[d]      = gridG->real->dx[d];
    gridL->fourierAl->dx[d] = gridG->fourierAl->dx[d];
    gridL->realAl->dx[d]    = gridG->realAl->dx[d];
  }

  // Also convenient to keep dealiased kperpSq in memory:
  gridL->fourier->xperpSq = mugy_alloc_real_ho(gridL->fourier->NxyTot);
  for (mint i=0; i<gridL->fourier->Nx[0]; i++) {
    for (mint j=0; j<gridL->fourier->Nx[1]; j++) {
      double kx = gridL->fourier->x[i];
      double ky = gridL->fourier->x[gridL->fourier->Nx[0]+j];
      gridL->fourier->xperpSq[i*gridL->fourier->Nx[1]+j] = kx*kx + ky*ky;
    }
  }

  // Set the type of each grid.
  gridL->real->type      = MUGY_REAL_GRID;
  gridL->fourier->type   = MUGY_FOURIER_GRID;
  gridL->realAl->type    = MUGY_REAL_GRID;
  gridL->fourierAl->type = MUGY_FOURIER_GRID;
}

void mugy_comms_terminate(struct mugy_comms *comms) {
  // Close MPI up, thoroughly.
  MPI_Barrier(comms->world->comm); // To avoid premature deallocations.

  struct mugy_comms_sub *scomm;  // Temp pointer for ease of notation.

  // Free world subcomm (DON'T FREE THE COMMUNICATOR).
  scomm = comms->world;
  free(scomm->decomp);  free(scomm->coord);
  free(scomm);

  // Free 4d subcomm.
  scomm = &comms->sub4d[0];
  free(scomm->decomp);  free(scomm->coord);
  MPI_Comm_free((MPI_Comm*) &scomm->comm);
  free(scomm);

  // Free 3d subcomm.
  scomm = &comms->sub3d[0];
  free(scomm->decomp);  free(scomm->coord);
  MPI_Comm_free((MPI_Comm*) &scomm->comm);
  free(comms->sub3d);

  // Free 2d subcomms.
  for (mint d=0; d<nDim; d++) {
    scomm = &comms->sub2d[d];
    free(scomm->decomp);  free(scomm->coord);
    MPI_Comm_free((MPI_Comm*) &scomm->comm);
  }
  free(comms->sub2d);

  // Free 1d subcomms.
  for (mint d=0; d<nDim+1; d++) {
    scomm = &comms->sub1d[d];
    free(scomm->decomp);  free(scomm->coord);
    MPI_Comm_free((MPI_Comm*) &scomm->comm);
  }
  free(comms->sub1d);

  free(comms);

  MPI_Finalize();
}
