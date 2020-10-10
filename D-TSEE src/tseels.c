#include<petsc.h>
PetscErrorCode tsee_error(Mat, Mat, Mat, PetscScalar, Mat, PetscInt, PetscInt, PetscScalar *);

PetscErrorCode tseels(Mat XXold, Mat XX, Mat Wp, Mat Wn, Mat ker, PetscScalar lambda, Mat P, PetscScalar ff, Mat G, PetscScalar *alpha, PetscInt cellcount, PetscInt d){

    PetscScalar rho=0.8, c=0.1, temp=0, e, *G_row, *P_row, tempsum;
    PetscInt i,j, istart, iend, *id;
    Mat tt;
    PetscErrorCode ierr;

    ierr = MatGetOwnershipRange(G, &istart, &iend); CHKERRQ(ierr);
    ierr = PetscMalloc1(d, &id); CHKERRQ(ierr);
    ierr = PetscMalloc1(d, &G_row); CHKERRQ(ierr);
    ierr = PetscMalloc1(d, &P_row); CHKERRQ(ierr);
    ierr = MatCreateSeqDense(MPI_COMM_SELF, cellcount, d, NULL, &tt); CHKERRQ(ierr);
    ierr = MatSeqDenseSetPreallocation(tt, NULL); CHKERRQ(ierr);

    for(i=0; i<d; i++){
        id[i]=i;
    }
    
    
    for(i=istart; i<iend; i++){
        ierr = MatGetValues(G, 1, &i, d, id, G_row); CHKERRQ(ierr);
        ierr = MatGetValues(P, 1, &i, d, id, P_row); CHKERRQ(ierr);
        for(j=0; j<d; j++){
            temp+=c*G_row[j]*P_row[j]; 
        }
    }
    ierr = MPI_Allreduce(&temp, &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(ierr);
    ierr = PetscFree(G_row); CHKERRQ(ierr);
    ierr = PetscFree(P_row); CHKERRQ(ierr);
    
    ierr = MatCopy(XXold, tt, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatAXPY(tt, *alpha, P, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = tsee_error(tt, Wp, Wn, lambda, ker, cellcount,d, &e); CHKERRQ(ierr);
    
    while(e > ff + *alpha * tempsum){
        *alpha = *alpha *rho;
        ierr = MatCopy(XX, tt, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = MatAXPY(tt, *alpha, P, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = tsee_error(tt, Wp, Wn, lambda, ker, cellcount, d, &e); CHKERRQ(ierr);
    }

    ierr = MatCopy(tt, XX, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = MatDestroy(&tt); CHKERRQ(ierr);
    ierr = PetscFree(id); CHKERRQ(ierr);

    return ierr;

}
