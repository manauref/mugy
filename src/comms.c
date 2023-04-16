/* mugy: mpi_comms

   Communication-related operations.
*/

#include "mh_comms.h"
#include <mpi.h>
#include "mh_utilities.h"
#include "mh_alloc.h"
#include "mh_io_tools.h"
#include <string.h>   // e.g. for memcpy.

struct mugy_comms *mugy_comms_init(mint argc, char *argv[]) {
  // Initialize MPI, get rank of this process and total number of processes.
  struct mugy_comms *comms = (struct mugy_comms *) malloc(sizeof(struct mugy_comms));
  comms->world.comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);  // Initialize MPI.
  MPI_Comm_rank(comms->world.comm, &comms->world.rank);  // Rank of this MPI process.
  MPI_Comm_size(comms->world.comm, &comms->world.size);  // Number of MPI processes.
  return comms;
}

void mugy_comms_sub_init(struct mugy_comms *comms, struct mugy_grid *grid, struct mugy_population *pop) {
  // Initialize the various sub-communicators needed.
  
  // Check the number of MPI processes is correct.
  if (prod_mint(grid->mpiProcs,nDim)*pop->mpiProcs != comms->world.size) {
    printf(" Number of MPI processes in input file (%d) differs from that in mpirun (%d).\n",
           prod_mint(grid->mpiProcs,nDim)*pop->mpiProcs, comms->world.size);
    abortSimulation(" Terminating...\n");
  }

  struct mugy_comms_sub *scomm;  // Temp pointer for ease of notation.

  // Finish initializing the world subcomm (coord initialized further down).
  scomm         = &comms->world; 
  scomm->dim    = nDim+1;
  scomm->coord  = alloc_mintArray_ho(scomm->dim);
  scomm->decomp = alloc_mintArray_ho(scomm->dim);
  for (mint d=0; d<nDim; d++) scomm->decomp[d] = grid->mpiProcs[d];
  scomm->decomp[nDim] = pop->mpiProcs;
  arrPrintS_mint(scomm->decomp,nDim+1, " MPI processes along X,Y,Z,s: ", "\n", comms->world.rank);
  r0printf("\n", comms->world.rank);

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
    scomm->decomp[d]      = comms->world.decomp[d];
    decompOrg[commOrg[d]] = comms->world.decomp[d];
    bcOrg[commOrg[d]]     = bc[d];
  }
  MPI_Cart_create(comms->world.comm, scomm->dim, decompOrg, bcOrg, reorder, (MPI_Comm*) &scomm->comm);
  MPI_Comm_size(scomm->comm, &scomm->size);
  MPI_Comm_rank(scomm->comm, &scomm->rank);
  MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  memcpy_mint(comms->world.coord, scomm->coord, nDim+1, host2host);  // Make coords the same in world.

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
      scomm->decomp[e] = comms->world.decomp[e];
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
      scomm->decomp[e] = comms->world.decomp[e];
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
    scomm->decomp[0] = comms->world.decomp[d];
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
  distribute1dDOFs(pop->mpiProcs, rank_s, pop->global.numSpecies, &pop->local.numSpecies, &pop->local.globalSpecOff);
  pop->local.globalMomOff = 0;
  for (mint s=0; s<pop->local.globalSpecOff; s++) pop->local.globalMomOff += pop->global.pars[s].numMoments;
  pop->local.pars = (struct mugy_species_pars*) calloc(pop->local.numSpecies, sizeof(struct mugy_species_pars));
  for (mint s=0; s<pop->local.numSpecies; s++) {
    pop->local.pars[s].numMoments = pop->global.pars[s+pop->local.globalSpecOff].numMoments;

    pop->local.pars[s].alpha      = mugy_alloc_real_ho(pop->local.pars[s].numMoments);
    pop->local.pars[s].nu         = mugy_alloc_real_ho(pop->local.pars[s].numMoments);
    pop->local.pars[s].hDiffOrder = mugy_alloc_real_ho(nDim);
    pop->local.pars[s].hDiff      = mugy_alloc_real_ho(nDim);
    pop->local.pars[s].kDiffMin   = mugy_alloc_real_ho(nDim);
    pop->local.pars[s].initAux    = mugy_alloc_real_ho(nDim);

    pop->local.pars[s].qCharge    = pop->global.pars[s+pop->local.globalSpecOff].qCharge   ;
    pop->local.pars[s].muMass     = pop->global.pars[s+pop->local.globalSpecOff].muMass    ;
    pop->local.pars[s].tau        = pop->global.pars[s+pop->local.globalSpecOff].tau       ;
    pop->local.pars[s].omSt       = pop->global.pars[s+pop->local.globalSpecOff].omSt      ;
    pop->local.pars[s].omd        = pop->global.pars[s+pop->local.globalSpecOff].omd       ;
    pop->local.pars[s].delta      = pop->global.pars[s+pop->local.globalSpecOff].delta     ;
    pop->local.pars[s].deltaPerp  = pop->global.pars[s+pop->local.globalSpecOff].deltaPerp ;
    pop->local.pars[s].eta        = pop->global.pars[s+pop->local.globalSpecOff].eta       ;
    memcpy(pop->local.pars[s].alpha, pop->global.pars[s+pop->local.globalSpecOff].alpha, pop->local.pars[s].numMoments*sizeof(real));
    memcpy(pop->local.pars[s].nu   , pop->global.pars[s+pop->local.globalSpecOff].nu   , pop->local.pars[s].numMoments*sizeof(real));
    pop->local.pars[s].delta0     = pop->global.pars[s+pop->local.globalSpecOff].delta0    ;
    memcpy(pop->local.pars[s].hDiffOrder, pop->global.pars[s+pop->local.globalSpecOff].hDiffOrder, nDim*sizeof(real));
    memcpy(pop->local.pars[s].hDiff     , pop->global.pars[s+pop->local.globalSpecOff].hDiff     , nDim*sizeof(real));
    memcpy(pop->local.pars[s].kDiffMin  , pop->global.pars[s+pop->local.globalSpecOff].kDiffMin  , nDim*sizeof(real));
    pop->local.pars[s].icOp       = pop->global.pars[s+pop->local.globalSpecOff].icOp      ;
    memcpy(pop->local.pars[s].initAux, pop->global.pars[s+pop->local.globalSpecOff].initAux, nDim*sizeof(real));
    pop->local.pars[s].initA      = pop->global.pars[s+pop->local.globalSpecOff].initA     ;
    pop->local.pars[s].noiseA     = pop->global.pars[s+pop->local.globalSpecOff].noiseA    ;
  }
  // Set the total number of moments.
  pop->local.numMomentsTot = 0;
  for (int s=0; s<pop->local.numSpecies; s++) pop->local.numMomentsTot += pop->local.pars[s].numMoments;

  // Distribute the real-space and Fourier-space points. 
  struct mugy_grid_ada *gridG = &grid->global;
  struct mugy_grid_ada *gridL = &grid->local;
  for (mint d=0; d<nDim; d++) {
    mint rank_i = comms->sub1d[d].rank;
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->deal.Nekx[d],
                     &gridL->deal.Nekx[d], &gridL->deal.globalOff[d]);
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->deal.dual.Nx[d],
                     &gridL->deal.dual.Nx[d], &gridL->deal.dual.globalOff[d]);
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->al.Nekx[d],
                     &gridL->al.Nekx[d], &gridL->al.globalOff[d]);
    distribute1dDOFs(grid->mpiProcs[d], rank_i, gridG->al.dual.Nx[d],
                     &gridL->al.dual.Nx[d], &gridL->al.dual.globalOff[d]);
  }
  gridL->deal.NekxTot      = prod_mint(gridL->deal.Nekx,nDim);
  gridL->deal.dual.NxTot   = prod_mint(gridL->deal.dual.Nx,nDim);
  gridL->deal.NekxyTot     = prod_mint(gridL->deal.Nekx,2);
  gridL->deal.dual.NxyTot  = prod_mint(gridL->deal.dual.Nx,2);
  gridL->al.NekxTot     = prod_mint(gridL->al.Nekx,nDim);
  gridL->al.dual.NxTot  = prod_mint(gridL->al.dual.Nx,nDim);
  gridL->al.NekxyTot    = prod_mint(gridL->al.Nekx,2);
  gridL->al.dual.NxyTot = prod_mint(gridL->al.dual.Nx,2);

  // Create local real-space and Fourier-space coordinate arrays.
  gridL->deal.kx     = mugy_alloc_real_ho(sum_mint(gridL->deal.Nekx, nDim));
  gridL->deal.dual.x = mugy_alloc_real_ho(sum_mint(gridL->deal.dual.Nx, nDim));
  gridL->al.kx     = mugy_alloc_real_ho(sum_mint(gridL->al.Nekx, nDim));
  gridL->al.dual.x = mugy_alloc_real_ho(sum_mint(gridL->al.dual.Nx, nDim));
  for (mint d=0; d<nDim; d++) {
    memcpy(getArray_real(gridL->deal.kx,gridL->deal.Nekx,d),
           getArray_real(gridG->deal.kx,gridG->deal.Nekx,d)+gridL->deal.globalOff[d], gridL->deal.Nekx[d]*sizeof(real));
    memcpy(getArray_real(gridL->deal.dual.x,gridL->deal.dual.Nx,d),
           getArray_real(gridG->deal.dual.x,gridG->deal.dual.Nx,d)+gridL->deal.dual.globalOff[d], gridL->deal.dual.Nx[d]*sizeof(real));
    memcpy(getArray_real(gridL->al.kx,gridL->al.Nekx,d),
           getArray_real(gridG->al.kx,gridG->al.Nekx,d)+gridL->al.globalOff[d], gridL->al.Nekx[d]*sizeof(real));
    memcpy(getArray_real(gridL->al.dual.x,gridL->al.dual.Nx,d),
           getArray_real(gridG->al.dual.x,gridG->al.dual.Nx,d)+gridL->al.dual.globalOff[d], gridL->al.dual.Nx[d]*sizeof(real));
  }

  // Copy global scalars into local grids.
  for (mint d=0; d<nDim; d++) {
    gridL->deal.kxMin[d]   = gridG->deal.kxMin[d];
    gridL->al.kxMin[d]     = gridG->al.kxMin[d];
    gridL->deal.dual.dx[d] = gridG->deal.dual.dx[d];
    gridL->al.dual.dx[d]   = gridG->al.dual.dx[d];
  }

  // Also convenient to keep dealiased kperpSq in memory:
  gridL->deal.kperpSq = mugy_alloc_real_ho(gridL->deal.NekxyTot);
  for (mint i=0; i<gridL->deal.Nekx[0]; i++) {
    for (mint j=0; j<gridL->deal.Nekx[1]; j++) {
      double kx = gridL->deal.kx[i];
      double ky = gridL->deal.kx[gridL->deal.Nekx[0]+j];
      gridL->deal.kperpSq[i*gridL->deal.Nekx[1]+j] = kx*kx + ky*ky;
    }
  }
}

void mugy_comms_terminate(struct mugy_comms *comms) {
  // Close MPI up, thoroughly.
  MPI_Barrier(comms->world.comm); // To avoid premature deallocations.

  struct mugy_comms_sub *scomm;  // Temp pointer for ease of notation.

  // Free world subcomm (DON'T FREE THE COMMUNICATOR).
  scomm = &comms->world;
  free(scomm->decomp);  free(scomm->coord);

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
