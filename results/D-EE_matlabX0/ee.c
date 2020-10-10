#include <petsc.h>
#include "DEE_myhead.h"

PetscErrorCode ee(Mat Wp, Mat Wn, PetscInt d, PetscScalar lambda, Mat XX, PetscInt allrowcount, PetscInt rowcount, PetscInt maxit, PetscScalar tol){

    PetscErrorCode ierr;
    PetscInt i, j, istart, iend, nlocal;
    PetscScalar alpha = 1.0;



    /* Set Lp4 and other variables required for iteration */
    Vec Dp;
    Mat Lp4;
    PetscScalar mDiagLp=0;

    ierr = VecCreateMPI(MPI_COMM_WORLD, rowcount, allrowcount, &Dp); CHKERRQ(ierr);
    ierr = MatGetRowSum(Wp, Dp); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Dp); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Dp); CHKERRQ(ierr);
    ierr = MatCreateDense(MPI_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, allrowcount, NULL, &Lp4); CHKERRQ(ierr);
    ierr = MatDuplicate(Wp, MAT_COPY_VALUES, &Lp4); CHKERRQ(ierr);

    ierr = MatScale(Lp4, -1.0); CHKERRQ(ierr);
    ierr = MatDiagonalSet(Lp4, Dp, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatScale(Lp4, 4.0); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Lp4, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Lp4, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    // define mDiagLp
    PetscInt ii=0;
    PetscScalar *array, tt_mDiagLp=0;
    ierr = VecGetSize(Dp, &nlocal); CHKERRQ(ierr);
    ierr = VecGetArray(Dp, &array); CHKERRQ(ierr);
    i = 0;
    while(i<nlocal){
        if(array[i]>0) break;
        i++;
    }
    ierr = MPI_Allreduce(&i, &ii, 1, MPI_INT,MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr); 
    ierr = MatGetOwnershipRange(Wp, &istart, &iend); CHKERRQ(ierr);  
    if(ii>=istart && ii<iend){
        tt_mDiagLp = 4*array[ii];
    }
    ierr = MPI_Allreduce(&tt_mDiagLp, &mDiagLp, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = VecRestoreArray(Dp, &array); CHKERRQ(ierr);
    

    // define Lp4_pd, which is an positive-definite Lp4 matrix and used to perform Cholesky factorization
    Mat Lp4_pd, Lp4_pd_elemental;
    ierr = MatCreateDense(MPI_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, allrowcount, NULL, &Lp4_pd); CHKERRQ(ierr);
    ierr = MatDuplicate(Lp4, MAT_COPY_VALUES, &Lp4_pd); CHKERRQ(ierr);
    ierr = MatShift(Lp4_pd, pow(10, -10)*mDiagLp); CHKERRQ(ierr);
    ierr = MatConvert(Lp4_pd, MATELEMENTAL, MAT_INITIAL_MATRIX, &Lp4_pd_elemental); CHKERRQ(ierr);
    ierr = MatSetOption(Lp4_pd_elemental, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatDestroy(&Lp4_pd); CHKERRQ(ierr);

    KSP ksp;
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, Lp4_pd_elemental, Lp4_pd_elemental); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    //define the initial error, ker
    PetscScalar e;
    Mat ker;
    ierr = MatCreateDense(MPI_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, allrowcount, NULL, &ker); CHKERRQ(ierr);
    ierr = ee_error(XX, Wp, Wn, lambda, ker, allrowcount, d, &e); CHKERRQ(ierr);


    /* set initial values of parameters and constant variables */
    PetscScalar difffronorm, oldfronorm;
    PetscInt convcrit = 1, itecount = 1;
    Mat G, P, L, XXold, XX_diff, WWn;
    Vec DDn;

    ierr = MatCreateDense(MPI_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, d, NULL, &G); CHKERRQ(ierr);
    ierr = MatCreateSeqDense(MPI_COMM_SELF, allrowcount, d, NULL, &P); CHKERRQ(ierr);
    ierr = MatCreateDense(MPI_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, allrowcount, NULL, &L); CHKERRQ(ierr);
    ierr = MatCreateSeqDense(MPI_COMM_SELF, allrowcount, d, NULL, &XXold); CHKERRQ(ierr);
    ierr = MatCreateSeqDense(MPI_COMM_SELF, allrowcount, d, NULL, &XX_diff); CHKERRQ(ierr);
    ierr = MatCreateDense(MPI_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, allrowcount, NULL, &WWn); CHKERRQ(ierr);
    ierr = VecCreateMPI(MPI_COMM_WORLD, rowcount, allrowcount, &DDn); 

    ierr = MatCopy(XX, XXold, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

    PetscScalar *Wnval, *kerval, *WWnval;
    ierr = PetscMalloc1(allrowcount, &WWnval); CHKERRQ(ierr);
    Vec X_temp, G_temp;
    PetscScalar *G_temp_array, *XX_temp_array, *XX_temp_col_array, *XX_temp_col_local_array;
    PetscInt vec_istart, vec_iend, vec_nlocal, *vec_id;

    ierr = VecCreateMPI(PETSC_COMM_WORLD, rowcount, allrowcount, &X_temp); CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, rowcount, allrowcount, &G_temp); CHKERRQ(ierr);
    ierr = PetscMalloc1(allrowcount, &XX_temp_col_array); CHKERRQ(ierr);
    ierr = PetscMalloc1(allrowcount, &vec_id); CHKERRQ(ierr);
    ierr = PetscMalloc1(allrowcount, &XX_temp_col_local_array); CHKERRQ(ierr);
    ierr = PetscMalloc1(allrowcount, &Wnval); CHKERRQ(ierr);
    ierr = PetscMalloc1(allrowcount, &kerval); CHKERRQ(ierr);


    PetscInt *id=NULL;
    PetscScalar *xx_val;
    ierr = PetscMalloc1(allrowcount, &id); CHKERRQ(ierr);
    ierr = PetscMalloc1(d, &xx_val); CHKERRQ(ierr);
    for(i=0;i<allrowcount;i++) id[i]=i;


    while(convcrit){
        
        for(i=istart; i<iend; i++){
            ierr = MatGetValues(Wn, 1, &i, allrowcount, id, Wnval); CHKERRQ(ierr);
            ierr = MatGetValues(ker, 1, &i, allrowcount, id, kerval); CHKERRQ(ierr);
            for(j=0; j<allrowcount; j++){
                WWnval[j] = lambda*Wnval[j]*kerval[j];
                if(j==i) WWnval[j] = 0;
            }
            ierr = MatSetValues(WWn, 1, &i, allrowcount, id, WWnval, INSERT_VALUES); CHKERRQ(ierr);   
        }
        ierr = MatAssemblyBegin(WWn, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(WWn, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatGetRowSum(WWn, DDn); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(DDn); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(DDn); CHKERRQ(ierr);

        ierr = MatScale(WWn, -1.0); CHKERRQ(ierr);
        ierr = MatDiagonalSet(WWn, DDn, ADD_VALUES); CHKERRQ(ierr);
        ierr = MatCopy(Lp4, L, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = MatAXPY(L, -4.0, WWn, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);        

        ierr = MatDenseGetArray(XX, &XX_temp_array); CHKERRQ(ierr);
        

        // compute gradient G=LX
        for(i=0; i<d; i++){
            for(j=0; j<allrowcount; j++){
                XX_temp_col_array[j] = XX_temp_array[i*allrowcount+j]; // XX_temp_col represents the i-th column of XX
            }
            // set values to MPI vector X_temp
            ierr = VecGetOwnershipRange(X_temp, &vec_istart, &vec_iend); CHKERRQ(ierr);
            ierr = VecGetLocalSize(X_temp, &vec_nlocal); CHKERRQ(ierr);
            for(j=0; j<vec_nlocal; j++){
                vec_id[j] = vec_istart+j;
                XX_temp_col_local_array[j] = XX_temp_col_array[vec_istart+j];
            }

            ierr = VecSetValues(X_temp, vec_nlocal, vec_id, XX_temp_col_local_array, INSERT_VALUES); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(X_temp); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(X_temp); CHKERRQ(ierr);
            // compute G_temp = L*X_temp
            ierr = MatMult(L, X_temp, G_temp); CHKERRQ(ierr);
            // access G_temp to g_temp_array
            ierr = VecGetArray(G_temp, &G_temp_array); CHKERRQ(ierr);
            // use G_temp_array to set values in i-th column of G
            ierr = MatSetValues(G, rowcount, vec_id, 1 ,&i, G_temp_array, INSERT_VALUES); CHKERRQ(ierr);
            ierr = VecRestoreArray(G_temp, &G_temp_array); CHKERRQ(ierr);
            
        }
    
        ierr = MatDenseRestoreArray(XX, &XX_temp_array); CHKERRQ(ierr);
        ierr = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


        // compute search direction
        ierr = SolveDir(ksp, G, P); CHKERRQ(ierr);

        // perform line search for steps
        ierr = eels(XXold, XX, Wp, Wn, ker, lambda, P, e, G, &alpha, allrowcount, d); CHKERRQ(ierr); CHKERRQ(ierr);
        ierr = ee_error(XX, Wp, Wn, lambda, ker, allrowcount, d, &e); CHKERRQ(ierr);


        // determine whether to stop
        ierr = MatCopy(XX, XX_diff, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = MatAXPY(XX_diff, -1.0, XXold, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

        ierr = MatNorm(XX_diff, NORM_FROBENIUS, &difffronorm); CHKERRQ(ierr);
        ierr = MatNorm(XXold, NORM_FROBENIUS, &oldfronorm); CHKERRQ(ierr);

        convcrit = (itecount<maxit) &&(difffronorm>tol*oldfronorm);
        
        // prepare for next iteration
        ierr = MatCopy(XX, XXold, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        itecount++;

    }
    ierr = VecDestroy(&X_temp); CHKERRQ(ierr);
    ierr = VecDestroy(&G_temp); CHKERRQ(ierr);
    ierr = PetscFree(XX_temp_col_array); CHKERRQ(ierr);
    ierr = PetscFree(vec_id); CHKERRQ(ierr);
    ierr = PetscFree(id); CHKERRQ(ierr);
    ierr = PetscFree(XX_temp_col_local_array); CHKERRQ(ierr);
    ierr = PetscFree(Wnval); CHKERRQ(ierr);
    ierr = PetscFree(kerval); CHKERRQ(ierr);
    ierr = MatDestroy(&Wp); CHKERRQ(ierr);
    ierr = MatDestroy(&Wn); CHKERRQ(ierr);
    ierr = MatDestroy(&ker); CHKERRQ(ierr);
    ierr = MatDestroy(&Lp4); CHKERRQ(ierr);
    ierr = MatDestroy(&Lp4_pd_elemental); CHKERRQ(ierr);
    ierr = MatDestroy(&G); CHKERRQ(ierr);
    ierr = MatDestroy(&P); CHKERRQ(ierr);
    ierr = MatDestroy(&L); CHKERRQ(ierr);
    ierr = MatDestroy(&XXold); CHKERRQ(ierr);
    ierr = MatDestroy(&XX_diff); CHKERRQ(ierr);
    ierr = MatDestroy(&WWn); CHKERRQ(ierr);
    ierr = VecDestroy(&Dp); CHKERRQ(ierr);
    ierr = VecDestroy(&DDn); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);


    
    return ierr;
}

