//
//  other_fun.c
//  parallel_test
//
//  Created by 安少坤 on 2019/6/22.
//  Copyright © 2019 Shaokun An. All rights reserved.
//
#include <stdlib.h>
#include <stdio.h>
#include <petsc.h>
PetscInt **int_get_mem(PetscInt rowcount, PetscInt colcount){
    PetscInt **data;
    data = (PetscInt**) malloc(rowcount*sizeof(PetscInt*));
    if(data==NULL){
        printf("Fail to allocate memeory for data!\n");
        exit(-1);
    }
    int i;
    for(i = 0; i < rowcount; i++){
        data[i] = (PetscInt*) malloc(colcount*sizeof(PetscInt));
        if(data[i]==NULL){
            printf("Fail to allocate memeory for data!\n");
            exit(-1);
        }
    }
    return data;
}

double **double_get_mem(PetscInt rowcount, PetscInt colcount){
    double **data;
    data = (double**) malloc(rowcount*sizeof(double*));
    if(data==NULL){
        printf("Fail to allocate memeory for data!\n");
        exit(-1);
    }
    int i;
    for(i = 0; i < rowcount; i++){
        data[i] = (double*) malloc(colcount*sizeof(double));
        if(data[i]==NULL){
            printf("Fail to allocate memeory for data!\n");
            exit(-1);
        }
    }
    return data;
}
