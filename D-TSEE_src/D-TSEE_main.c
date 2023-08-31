// 
//  main.c
//  parallelized Elastic Embedding
//
//  Created by 安少坤 on 2019/3/25.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <petsc.h>
#include "myhead.h"
#include <math.h>


int main(int argc, char **args){

    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    PetscInt i,j;


    ierr = PetscInitialize(&argc, &args, (char*)0, (char*)0); if(ierr) return ierr;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    char filename1[100];
    char timefilename[100];
    PetscBool flg1, flg2;
    ierr = PetscOptionsGetString(NULL, NULL, "-inputfile", filename1, sizeof(filename1), &flg1); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL, NULL, "-timefile", timefilename, sizeof(timefilename), &flg2); CHKERRQ(ierr);
    if(flg1 & flg2){
        FILE *fp;
        fp = fopen(filename1, "r");
        if(fp==NULL){
            ierr = PetscPrintf(PETSC_COMM_WORLD, "Error: Fail to read Data!\n"); CHKERRQ(ierr);
            fclose(fp);
            return ierr;
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Read data done!\n"); CHKERRQ(ierr);
        PetscInt colcount, allrowcount; //allrowcount is the number of cells, colcount is the number of features
        ierr = GetCsvColCount(fp, &colcount); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "feature count = %d\n", colcount); CHKERRQ(ierr);
        ierr = GetCsvRowCount(fp, &allrowcount); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "cell count = %d\n", allrowcount); CHKERRQ(ierr);
        
        
        // read original data in all processes
        double **data=NULL;
        data = double_get_mem(allrowcount, colcount);
        ierr = ReadCsvData(data, fp); CHKERRQ(ierr);
        ierr = NormalizeData(data, allrowcount, colcount); CHKERRQ(ierr);

        FILE *fptime;
        fptime = fopen(timefilename, "r");
        if(fptime==NULL){
            ierr = PetscPrintf(PETSC_COMM_WORLD, "Error: Fail to read time file!\n"); CHKERRQ(ierr);
            fclose(fptime);
            return ierr;
        }

        double **time=NULL;
        time = double_get_mem(allrowcount, 1);
        ierr = ReadCsvData(time, fptime); CHKERRQ(ierr);
        ierr = NormalizeData(time, allrowcount, 1); CHKERRQ(ierr);

        PetscInt rowcount;//rowcount is the number of cells distributed in each process

        // compute Wp in each process 

        Mat Wp, Wn;
        ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, allrowcount, allrowcount, NULL, &Wp); CHKERRQ(ierr);
        ierr = MatGetLocalSize(Wp, &rowcount, NULL); CHKERRQ(ierr);
        ierr = MatCreateDense(PETSC_COMM_WORLD, rowcount, PETSC_DECIDE, allrowcount, allrowcount, NULL, &Wn); CHKERRQ(ierr);

        ierr = ea(data, 20, allrowcount, rowcount, colcount, Wp); CHKERRQ(ierr);

        // cpmpute Wn in each process

        PetscScalar *temp_Wn, *temp_time_Wn;
        PetscInt istart, iend, *id;
        PetscScalar beta = 10.0;
        ierr = PetscOptionsGetReal(NULL, NULL, "-beta", &beta, NULL); CHKERRQ(ierr);
        ierr = PetscMalloc1(allrowcount, &temp_Wn); CHKERRQ(ierr);
        ierr = PetscMalloc1(allrowcount, &temp_time_Wn); CHKERRQ(ierr);
        ierr = PetscMalloc1(allrowcount, &id); CHKERRQ(ierr);
        ierr = MatGetOwnershipRange(Wn, &istart, &iend); CHKERRQ(ierr);
        for(i=0;i<allrowcount;i++) id[i]=i;
        for(i=istart;i<iend; i++){
            for(j=0; j<allrowcount; j++){
                ierr = sqdist(data[i], data[j], colcount, temp_Wn+j); CHKERRQ(ierr);
                ierr = absdist(time[i],time[j], 1, temp_time_Wn+j); CHKERRQ(ierr);
                temp_Wn[j] += beta*temp_time_Wn[j];
            }
            ierr = MatSetValues(Wn, 1, &i, allrowcount, id, temp_Wn, INSERT_VALUES); CHKERRQ(ierr);
        }
        ierr = MatAssemblyBegin(Wn, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Wn, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = PetscFree(temp_Wn); CHKERRQ(ierr);
        ierr = PetscFree(id); CHKERRQ(ierr);
        ierr = PetscFree(temp_time_Wn); CHKERRQ(ierr);


        // make Wp and Wn to be symmetric and zero diagonal entries
        Vec zerovec;
        ierr = VecCreateMPI(PETSC_COMM_WORLD, rowcount, allrowcount, &zerovec); CHKERRQ(ierr);
        ierr = VecSet(zerovec, 0.0); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(zerovec); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(zerovec); CHKERRQ(ierr);
        // zero diagonal entries
        ierr = MatDiagonalSet(Wp, zerovec, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatDiagonalSet(Wn, zerovec, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecDestroy(&zerovec); CHKERRQ(ierr);
        // symmetric mat
        Mat Wt;
        ierr = MatConvert(Wp, MATSAME, MAT_INITIAL_MATRIX, &Wt); CHKERRQ(ierr);
        ierr = MatTranspose(Wp, MAT_REUSE_MATRIX, &Wt); CHKERRQ(ierr);
        ierr = MatAXPY(Wp, 1.0, Wt, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = MatScale(Wp, 0.5); CHKERRQ(ierr);
        ierr = MatTranspose(Wn, MAT_REUSE_MATRIX, &Wt); CHKERRQ(ierr);
        ierr = MatAXPY(Wn, 1.0, Wt, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = MatScale(Wn, 0.5); CHKERRQ(ierr);
        ierr = MatDestroy(&Wt); CHKERRQ(ierr);


        // normalize entries by dividing the summ of all entries for each entry
        Vec Wp_sum, Wn_sum;
        PetscScalar wp_sum_sum, wn_sum_sum;
        ierr = VecCreateMPI(PETSC_COMM_WORLD, rowcount, allrowcount, &Wp_sum); CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, rowcount, allrowcount, &Wn_sum); CHKERRQ(ierr);
        ierr = MatGetRowSum(Wp,Wp_sum); CHKERRQ(ierr);
        ierr = MatGetRowSum(Wn, Wn_sum); CHKERRQ(ierr);
        ierr = VecSum(Wp_sum, &wp_sum_sum); CHKERRQ(ierr);
        ierr = VecSum(Wn_sum, &wn_sum_sum); CHKERRQ(ierr);
        ierr = MatScale(Wp, 1/(wp_sum_sum/pow(10,6))); CHKERRQ(ierr);
        ierr = MatScale(Wn, 1/(wn_sum_sum/pow(10,6))); CHKERRQ(ierr);
        ierr = VecDestroy(&Wp_sum); CHKERRQ(ierr);
        ierr = VecDestroy(&Wn_sum); CHKERRQ(ierr);


        // compute EE in each processor
        PetscInt d=2, maxit=100;
        PetscScalar lambda=10.0, tol = pow(10,-3);
        ierr = PetscOptionsGetReal(NULL, NULL, "-lambda", &lambda, NULL); CHKERRQ(ierr);

        double** XX_val;
        XX_val = (double**)malloc(allrowcount*sizeof(double*));
        for(i=0; i<allrowcount; i++){
            XX_val[i] = (double*)malloc(d*sizeof(double));
        }

        double x0_mean=0, x0_sd=0;
        for(i=0; i<allrowcount; i++){
            for(j=0; j<d; j++){
                XX_val[i][j] = rand();
                x0_mean += XX_val[i][j];
            }
        }
        x0_mean = x0_mean/(allrowcount*d);
        for(i=0; i<allrowcount; i++){
            for(j=0; j<d; j++){
                x0_sd += (XX_val[i][j]-x0_mean)*(XX_val[i][j]-x0_mean);
            }
        }
        x0_sd = x0_sd/(allrowcount*d);
        x0_sd = sqrt(x0_sd);
        for(i=0; i<allrowcount; i++){
            for(j=0; j<d; j++){
                XX_val[i][j] = (XX_val[i][j]-x0_mean)/x0_sd*1e-5;
            }
        }


        Mat XX;
        ierr = MatCreateSeqDense(PETSC_COMM_SELF, allrowcount, d, NULL, &XX); CHKERRQ(ierr);
        PetscInt *idx;
        ierr = PetscMalloc1(d, &idx); CHKERRQ(ierr);
        for(i=0; i<d; i++) idx[i] = i;
        for(i=0; i<allrowcount; i++){
            ierr = MatSetValues(XX, 1, &i, 2, idx, XX_val[i], INSERT_VALUES); CHKERRQ(ierr);
        }
        ierr = MatAssemblyBegin(XX, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(XX, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = PetscFree(idx); CHKERRQ(ierr);
        free(XX_val);

        ierr = tsee(Wp, Wn, d, lambda, XX, allrowcount, rowcount, maxit, tol); CHKERRQ(ierr);


        // write data into csv file

        if(rank==0){
            double *xx;
            PetscScalar tt;
            xx = (double*)malloc(allrowcount*d*sizeof(double));
            for(i=0; i<allrowcount; i++){
                for(j=0; j<d; j++){
                    ierr = MatGetValues(XX, 1, &i, 1, &j, &tt); CHKERRQ(ierr);
                    xx[i*d+j] = tt;
                }
            }    
            char outname[100];
            PetscBool flg3;
            ierr = PetscOptionsGetString(NULL, NULL, "-outfile", outname, sizeof(outname), &flg3); CHKERRQ(ierr);
            if(flg3){
                ierr = WriteCsvData(outname, xx, allrowcount, d); CHKERRQ(ierr);
            }else{
                ierr = PetscPrintf(PETSC_COMM_WORLD, "Error: Please input the output file name with -outfile\n"); CHKERRQ(ierr);
            }

            free(xx);
        }
    

        ierr = MatDestroy(&XX); CHKERRQ(ierr); 

    }else{
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Error: Not found the input data file!\n"); CHKERRQ(ierr);
    }
    
    
    ierr = PetscFinalize();
    
    return ierr;

}
