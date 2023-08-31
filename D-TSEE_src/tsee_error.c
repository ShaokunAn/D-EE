#include <petsc.h>
PetscErrorCode tsee_error(Mat XX, Mat Wp, Mat Wn, PetscScalar lambda, Mat ker, PetscInt cellcount, PetscInt d, PetscScalar* e){

    PetscErrorCode ierr;
    PetscInt istart, iend, i, j, k, *id;
    PetscScalar *XX_val, sqd_val=0, *kerval, tt;

    ierr = MatGetOwnershipRange(ker, &istart, &iend); CHKERRQ(ierr);
    ierr = MatDenseGetArray(XX, &XX_val); CHKERRQ(ierr);
    ierr = PetscMalloc1(cellcount, &id); CHKERRQ(ierr);
    
    for(i=0; i<cellcount; i++) 
        id[i] = i;

    // compute ker first
    ierr = PetscMalloc1(cellcount, &kerval); CHKERRQ(ierr);
    for(i=istart; i<iend; i++){
        for(j=0; j<cellcount; j++){
            sqd_val = 0;
            for(k=0; k<d; k++){
                tt = XX_val[k*cellcount+j]-XX_val[k*cellcount+i];
                sqd_val = sqd_val + pow(tt,2);
            }
            kerval[j] = exp(-sqd_val); 
        }
        ierr = MatSetValues(ker, 1, &i, cellcount, id, kerval, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(ker, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(ker, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = PetscFree(kerval); CHKERRQ(ierr);
    

    // compute e value 
    PetscScalar e_temp = 0;
    PetscScalar *Wp_rowval, *Wn_rowval;
    ierr = PetscMalloc1(cellcount, &Wp_rowval); CHKERRQ(ierr);
    ierr = PetscMalloc1(cellcount, &Wn_rowval); CHKERRQ(ierr);
    
    for(i=istart; i<iend; i++){
        ierr = MatGetValues(Wp, 1, &i, cellcount, id, Wp_rowval); CHKERRQ(ierr);
        ierr = MatGetValues(Wn, 1, &i, cellcount, id, Wn_rowval); CHKERRQ(ierr);
        for(j=i+1; j<cellcount; j++){
            sqd_val = 0;
            for(k=0; k<d; k++){
                tt = XX_val[k*cellcount+j]-XX_val[k*cellcount+i];
                sqd_val = sqd_val + pow(tt,2);
            }
            e_temp+=2*(Wp_rowval[j]*sqd_val+lambda*Wn_rowval[j]*exp(-sqd_val));
        }
    }

    ierr = MPI_Allreduce(&e_temp,e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); CHKERRQ(ierr);

    
    // free memory
    ierr = PetscFree(id); CHKERRQ(ierr);
    ierr = PetscFree(Wp_rowval); CHKERRQ(ierr);
    ierr = PetscFree(Wn_rowval); CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(XX, &XX_val); CHKERRQ(ierr);
    return ierr;    

}
