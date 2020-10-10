//
//  other_fun.c
//  parallel_test
//
//  Created by 安少坤 on 2019/6/22.
//  Copyright © 2019 Shaokun An. All rights reserved.
//
#include <stdlib.h>
#include <stdio.h>
int **int_get_mem(int rowcount, int colcount){
    int **data;
    data = (int**) malloc(rowcount*sizeof(int*));
    if(data==NULL){
        printf("Fail to allocate memeory for data!\n");
        exit(-1);
    }
    for(int i = 0; i < rowcount; i++){
        data[i] = (int*) malloc(colcount*sizeof(int));
        if(data[i]==NULL){
            printf("Fail to allocate memeory for data!\n");
            exit(-1);
        }
    }
    return data;
}

double **double_get_mem(int rowcount, int colcount){
    double **data;
    data = (double**) malloc(rowcount*sizeof(double*));
    if(data==NULL){
        printf("Fail to allocate memeory for data!\n");
        exit(-1);
    }
    for(int i = 0; i < rowcount; i++){
        data[i] = (double*) malloc(colcount*sizeof(double));
        if(data[i]==NULL){
            printf("Fail to allocate memeory for data!\n");
            exit(-1);
        }
    }
    return data;
}
