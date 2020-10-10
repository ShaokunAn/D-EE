//
//  NormalizeData.c
//  parallel_test
//
//  Created by 安少坤 on 2019/3/26.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>

#define MAX_LINE_SIZE 1000000

PetscErrorCode FindMin(double *num, PetscInt len, double *min);
PetscErrorCode FindMax(double *num, PetscInt len, double *max);


PetscErrorCode NormalizeData(double** data, PetscInt rowcount, PetscInt colcount){
    
    PetscErrorCode ierr;
    double min, max;
    double* RowMin = (double*) malloc(rowcount*sizeof(double));
    double* RowMax = (double*) malloc(rowcount*sizeof(double));

    int i, j;
    for (i = 0; i<rowcount; i++) {
        ierr = FindMin(*(data+i), colcount, RowMin+i); CHKERRQ(ierr);
        ierr = FindMax(*(data+i), colcount, RowMax+i); CHKERRQ(ierr);
    }
    ierr = FindMin(RowMin, rowcount, &min); CHKERRQ(ierr);
    ierr = FindMax(RowMax, rowcount, &max);
//    printf("max value of data is %f\n", max);
    free(RowMin);
    free(RowMax);
    
    double denomintor = max-min;
    for (i = 0; i<rowcount; i++) {
        for (j = 0; j<colcount; j++) {
            data[i][j] = (data[i][j]-min)/denomintor;
        }
    }
    return ierr;

}

PetscErrorCode FindMin(double *num, PetscInt len, double *min){
    PetscErrorCode ierr = 0;
    *min=*num;
    int i;
    for(i=0;i<len;i++){
        if(*min>*(num+i)) *min = *(num+i);
    }
    return ierr;
}

PetscErrorCode FindMax(double *num, PetscInt len, double *max){
    PetscErrorCode ierr = 0;
    *max = *num;
    int i;
    for(i=0; i<len;i++){
        if(*max<*(num+i)) *max = *(num+i);
    }
    return ierr;
}

