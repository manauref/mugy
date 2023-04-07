/* mugy: mpi_comms

   Communication-related operations.
*/

#include "mh_comms.h"
#include <string.h>   // e.g. for memcpy.

void comms_init(struct mugy_comms *comms, mint argc, char *argv[]) {
  // Initialize MPI, get rank of this process and total number of processes.
  comms->world.comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);  // Initialize MPI.
  MPI_Comm_rank(comms->world.comm, &comms->world.rank);  // Rank of this MPI process.
  MPI_Comm_size(comms->world.comm, &comms->world.size);  // Number of MPI processes.
}

void comms_sub_init(struct mugy_comms *comms, struct mugy_grid grid, struct mugy_population pop) {
  // Initialize the various sub-communicators needed.
  
  // Check the number of MPI processes is correct.
  if (prod_mint(grid.mpiProcs,nDim)*pop.mpiProcs != comms->world.size) {
    printf(" Number of MPI processes in input file (%d) differs from that in mpirun (%d).\n",
           prod_mint(grid.mpiProcs,nDim)*pop.mpiProcs, comms->world.size);
    abortSimulation(" Terminating...\n");
  }

  struct mugy_comms_sub *scomm;  // Temp pointer for ease of notation.

  // Finish initializing the world subcomm (coord initialized further down).
  scomm         = &comms->world; 
  scomm->dim    = nDim+1;
  scomm->coord  = alloc_mintArray_ho(scomm->dim);
  scomm->decomp = alloc_mintArray_ho(scomm->dim);
  for (mint d=0; d<nDim; d++) scomm->decomp[d] = grid.mpiProcs[d];
  scomm->decomp[nDim] = pop.mpiProcs;
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
  MPI_Cart_create(comms->world.comm, scomm->dim, decompOrg, bcOrg, reorder, &scomm->comm);
  MPI_Comm_size(scomm->comm, &scomm->size);
  MPI_Comm_rank(scomm->comm, &scomm->rank);
  MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  memcpy_mint(comms->world.coord, scomm->coord, nDim+1, host2host);  // Make coords the same in world.

  // 3D subcommunicator (xyz space, no need for xys, xzs, yzs comms).
  comms->sub3d = (struct mugy_comms_sub *) malloc(sizeof(struct mugy_comms_sub));
  for (mint d=nDim; d<nDim+1; d++) {
    mint remain[nDim+1] = {true,true,true,true};
    remain[commOrg[d]] = false;

    scomm            = &comms->sub3d[0]; 
    scomm->dim       = 3;
    scomm->coord     = alloc_mintArray_ho(scomm->dim);
    scomm->decomp    = alloc_mintArray_ho(scomm->dim);
    for (mint e=0; e<scomm->dim; e++)
      scomm->decomp[e] = comms->world.decomp[e];
    MPI_Cart_sub(comms->sub4d[0].comm, remain, &scomm->comm);
    MPI_Comm_size(scomm->comm, &scomm->size);
    MPI_Comm_rank(scomm->comm, &scomm->rank);
    MPI_Cart_coords(scomm->comm, scomm->rank, scomm->dim, scomm->coord);
  }

  // 2D subcommunicators (no need for xs,ys,zs subcomms at this point).
  comms->sub2d = (struct mugy_comms_sub *) malloc(nDim*sizeof(struct mugy_comms_sub));
  for (mint d=nDim-1; d>-1; d--) {
    mint remain[nDim+1] = {true,true,true,true};
    remain[commOrg[d]] = false;  remain[commOrg[nDim]] = false;

    scomm            = &comms->sub2d[d]; 
    scomm->dim       = 2;
    scomm->coord     = alloc_mintArray_ho(scomm->dim);
    scomm->decomp    = alloc_mintArray_ho(scomm->dim);
    for (mint e=0; e<scomm->dim; e++)
      scomm->decomp[e] = comms->world.decomp[e];
    MPI_Cart_sub(comms->sub4d[0].comm, remain, &scomm->comm);
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
    MPI_Cart_sub(comms->sub4d[0].comm, remain, &scomm->comm);
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

void distributeDOFs(struct mugy_comms comms, struct mugy_grid globalGrid, struct mugy_population globalPop,
  struct mugy_grid *localGrid, struct mugy_population *localPop) {
  // Distribute s,Z,X,Y amongst MPI processes.
  
  // Distribute the species.
  mint rank_s = comms.sub1d[nDim].rank;
  distribute1dDOFs(globalPop.mpiProcs, rank_s, globalPop.numSpecies, &localPop->numSpecies, &localPop->globalSpecOff);
  localPop->globalMomOff = 0;
  for (mint s=0; s<localPop->globalSpecOff; s++) localPop->globalMomOff += globalPop.spec[s].numMoments;
  localPop->spec = (struct mugy_species*) calloc(localPop->numSpecies, sizeof(struct mugy_species));
  for (mint s=0; s<localPop->numSpecies; s++) {
    localPop->spec[s].numMoments = globalPop.spec[s+localPop->globalSpecOff].numMoments;

    localPop->spec[s].alpha      = alloc_realArray_ho(localPop->spec[s].numMoments);
    localPop->spec[s].nu         = alloc_realArray_ho(localPop->spec[s].numMoments);
    localPop->spec[s].hDiffOrder = alloc_realArray_ho(nDim);
    localPop->spec[s].hDiff      = alloc_realArray_ho(nDim);
    localPop->spec[s].kDiffMin   = alloc_realArray_ho(nDim);
    localPop->spec[s].initAux    = alloc_realArray_ho(nDim);

    localPop->spec[s].qCharge    = globalPop.spec[s+localPop->globalSpecOff].qCharge   ;
    localPop->spec[s].muMass     = globalPop.spec[s+localPop->globalSpecOff].muMass    ;
    localPop->spec[s].tau        = globalPop.spec[s+localPop->globalSpecOff].tau       ;
    localPop->spec[s].omSt       = globalPop.spec[s+localPop->globalSpecOff].omSt      ;
    localPop->spec[s].omd        = globalPop.spec[s+localPop->globalSpecOff].omd       ;
    localPop->spec[s].delta      = globalPop.spec[s+localPop->globalSpecOff].delta     ;
    localPop->spec[s].deltaPerp  = globalPop.spec[s+localPop->globalSpecOff].deltaPerp ;
    localPop->spec[s].eta        = globalPop.spec[s+localPop->globalSpecOff].eta       ;
    memcpy(localPop->spec[s].alpha, globalPop.spec[s+localPop->globalSpecOff].alpha, localPop->spec[s].numMoments*sizeof(real));
    memcpy(localPop->spec[s].nu   , globalPop.spec[s+localPop->globalSpecOff].nu   , localPop->spec[s].numMoments*sizeof(real));
    localPop->spec[s].delta0     = globalPop.spec[s+localPop->globalSpecOff].delta0    ;
    memcpy(localPop->spec[s].hDiffOrder, globalPop.spec[s+localPop->globalSpecOff].hDiffOrder, nDim*sizeof(real));
    memcpy(localPop->spec[s].hDiff     , globalPop.spec[s+localPop->globalSpecOff].hDiff     , nDim*sizeof(real));
    memcpy(localPop->spec[s].kDiffMin  , globalPop.spec[s+localPop->globalSpecOff].kDiffMin  , nDim*sizeof(real));
    localPop->spec[s].icOp       = globalPop.spec[s+localPop->globalSpecOff].icOp      ;
    memcpy(localPop->spec[s].initAux, globalPop.spec[s+localPop->globalSpecOff].initAux, nDim*sizeof(real));
    localPop->spec[s].initA      = globalPop.spec[s+localPop->globalSpecOff].initA     ;
    localPop->spec[s].noiseA     = globalPop.spec[s+localPop->globalSpecOff].noiseA    ;
  }
  // Set the total number of moments.
  localPop->numMomentsTot = 0;
  for (int s=0; s<localPop->numSpecies; s++) localPop->numMomentsTot += localPop->spec[s].numMoments;

  // Distribute the real-space and Fourier-space points. 
  for (mint d=0; d<nDim; d++) {
    mint rank_i = comms.sub1d[d].rank;
    distribute1dDOFs(globalGrid.mpiProcs[d], rank_i, globalGrid.fG.Nekx[d],
                     &localGrid->fG.Nekx[d], &localGrid->fG.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], rank_i, globalGrid.fG.dual.Nx[d],
                     &localGrid->fG.dual.Nx[d], &localGrid->fG.dual.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], rank_i, globalGrid.fGa.Nekx[d],
                     &localGrid->fGa.Nekx[d], &localGrid->fGa.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], rank_i, globalGrid.fGa.dual.Nx[d],
                     &localGrid->fGa.dual.Nx[d], &localGrid->fGa.dual.globalOff[d]);
  }
  localGrid->fG.NekxTot      = prod_mint(localGrid->fG.Nekx,nDim);
  localGrid->fG.dual.NxTot   = prod_mint(localGrid->fG.dual.Nx,nDim);
  localGrid->fG.NekxyTot     = prod_mint(localGrid->fG.Nekx,2);
  localGrid->fG.dual.NxyTot  = prod_mint(localGrid->fG.dual.Nx,2);
  localGrid->fGa.NekxTot     = prod_mint(localGrid->fGa.Nekx,nDim);
  localGrid->fGa.dual.NxTot  = prod_mint(localGrid->fGa.dual.Nx,nDim);
  localGrid->fGa.NekxyTot    = prod_mint(localGrid->fGa.Nekx,2);
  localGrid->fGa.dual.NxyTot = prod_mint(localGrid->fGa.dual.Nx,2);

  // Create local real-space and Fourier-space coordinate arrays.
  localGrid->fG.kx     = alloc_realArray_ho(sum_mint(localGrid->fG.Nekx, nDim));
  localGrid->fG.dual.x = alloc_realArray_ho(sum_mint(localGrid->fG.dual.Nx, nDim));
  localGrid->fGa.kx     = alloc_realArray_ho(sum_mint(localGrid->fGa.Nekx, nDim));
  localGrid->fGa.dual.x = alloc_realArray_ho(sum_mint(localGrid->fGa.dual.Nx, nDim));
  for (mint d=0; d<nDim; d++) {
    memcpy(getArray_real(localGrid->fG.kx,localGrid->fG.Nekx,d),
           getArray_real(globalGrid.fG.kx,globalGrid.fG.Nekx,d)+localGrid->fG.globalOff[d], localGrid->fG.Nekx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fG.dual.x,localGrid->fG.dual.Nx,d),
           getArray_real(globalGrid.fG.dual.x,globalGrid.fG.dual.Nx,d)+localGrid->fG.dual.globalOff[d], localGrid->fG.dual.Nx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fGa.kx,localGrid->fGa.Nekx,d),
           getArray_real(globalGrid.fGa.kx,globalGrid.fGa.Nekx,d)+localGrid->fGa.globalOff[d], localGrid->fGa.Nekx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fGa.dual.x,localGrid->fGa.dual.Nx,d),
           getArray_real(globalGrid.fGa.dual.x,globalGrid.fGa.dual.Nx,d)+localGrid->fGa.dual.globalOff[d], localGrid->fGa.dual.Nx[d]*sizeof(real));
  }

  // Copy global scalars into local grids.
  for (mint d=0; d<nDim; d++) {
    localGrid->fG.kxMin[d]    = globalGrid.fG.kxMin[d];
    localGrid->fGa.kxMin[d]   = globalGrid.fGa.kxMin[d];
    localGrid->fG.dual.dx[d]  = globalGrid.fG.dual.dx[d];
    localGrid->fGa.dual.dx[d] = globalGrid.fGa.dual.dx[d];
  }

}

void comms_terminate(struct mugy_comms *comms) {
  // Close MPI up, thoroughly.
  MPI_Barrier(comms->world.comm); // To avoid premature deallocations.

  struct mugy_comms_sub *scomm;  // Temp pointer for ease of notation.

  // Free world subcomm (DON'T FREE THE COMMUNICATOR).
  scomm = &comms->world;
  free(scomm->decomp);  free(scomm->coord);

  // Free 4d subcomm.
  scomm = &comms->sub4d[0];
  free(scomm->decomp);  free(scomm->coord);
  MPI_Comm_free(&scomm->comm);
  free(scomm);

  // Free 3d subcomm.
  scomm = &comms->sub3d[0];
  free(scomm->decomp);  free(scomm->coord);
  MPI_Comm_free(&scomm->comm);
  free(comms->sub3d);

  // Free 2d subcomms.
  for (mint d=0; d<nDim; d++) {
    scomm = &comms->sub2d[d];
    free(scomm->decomp);  free(scomm->coord);
    MPI_Comm_free(&scomm->comm);
  }
  free(comms->sub2d);

  // Free 1d subcomms.
  for (mint d=0; d<nDim+1; d++) {
    scomm = &comms->sub1d[d];
    free(scomm->decomp);  free(scomm->coord);
    MPI_Comm_free(&scomm->comm);
  }
  free(comms->sub1d);

  MPI_Finalize();
}
