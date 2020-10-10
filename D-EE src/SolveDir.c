#include<petsc.h>
#include<stdlib.h>

PetscErrorCode SolveDir(KSP ksp, Mat G, Mat P){
    /* G is Mat, need to be transformed to vec to use KSP, and the same as P*/

    PetscErrorCode ierr;
    PetscInt N, d, i, j, *id, istart, iend, nlocal;
    int *recvcounts, *displs;
    Vec x_P, rhs_G;
    PetscScalar *val, *all_val;
    Mat negG;
    PetscMPIInt rank, size;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    ierr = MatConvert(G, MATDENSE, MAT_INITIAL_MATRIX, &negG); CHKERRQ(ierr);
    ierr = MatScale(negG, -1.0); CHKERRQ(ierr);

    ierr = MatGetSize(negG, &N, &d); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(negG, &istart, &iend); CHKERRQ(ierr);
    ierr = MatGetLocalSize(negG, &nlocal, NULL); CHKERRQ(ierr);
    

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N, &x_P); CHKERRQ(ierr);
    ierr = VecDuplicate(x_P, &rhs_G); CHKERRQ(ierr);

    ierr = PetscMalloc1(N, &id); CHKERRQ(ierr);
    ierr = PetscMalloc1(N, &all_val); CHKERRQ(ierr);
    recvcounts = (int*) malloc(size*sizeof(int));
    displs = (int*) malloc(size*sizeof(int));

    for(i=0; i<N; i++){
        id[i]=i;
    }

    ierr = MPI_Allgather(&nlocal, 1, MPI_INT, recvcounts, 1, MPI_INT,MPI_COMM_WORLD); CHKERRQ(ierr);
    displs[0] = 0;
    int dispi;
    for(i=1; i<size; i++){
        dispi = 0;
        for(j=0; j<i; j++){
            dispi+=recvcounts[j];
        }
        displs[i] = dispi;
    }



    // for(i=0; i<size; i++){
    //     displs[i] = i*N/size;
    //     if(rank==size-1){
    //         recvcounts[i] = N/size+N%size;
    //     }else{
    //         recvcounts[i] = N/size;
    //     }
    // }

    /* 
        turn G and P into corresponding vectors, respectively.
        and solve a sequence of linear system with different RHS 
    */

    for(i=0; i<d; i++){
        ierr = MatGetColumnVector(negG, rhs_G, i); CHKERRQ(ierr); // this routine may be SLOW!!!
        ierr = KSPSolve(ksp, rhs_G, x_P); CHKERRQ(ierr);
        ierr = VecGetArray(x_P, &val); CHKERRQ(ierr);
        ierr = MPI_Allgatherv(val, nlocal, MPI_DOUBLE, all_val, recvcounts, displs, MPI_DOUBLE, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MatSetValues(P, N, id, 1, &i, all_val, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecRestoreArray(x_P, &val); CHKERRQ(ierr);
    }

    /* free memory*/
    
    ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    free(recvcounts);
    free(displs);
    ierr = PetscFree(id); CHKERRQ(ierr);
    ierr = PetscFree(all_val); CHKERRQ(ierr);
    ierr = VecDestroy(&x_P); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs_G); CHKERRQ(ierr);
    ierr = MatDestroy(&negG); CHKERRQ(ierr);

    return ierr;

}
