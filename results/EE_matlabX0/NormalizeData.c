//
//  NormalizeData.c
//  parallel_test
//
//  Created by 安少坤 on 2019/3/26.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE_SIZE 1000000

double FindMin(double* num, int len);
double FindMax(double* num, int len);


void NormalizeData(double** data, int rowcount, int colcount){
    
    double min, max;
    double* RowMin = (double*) malloc(rowcount*sizeof(double));
    double* RowMax = (double*) malloc(rowcount*sizeof(double));
    
    for (int i = 0; i<rowcount; i++) {
        *(RowMin+i) = FindMin(*(data+i), colcount);
        *(RowMax+i) = FindMax(*(data+i), colcount);
    }
    min = FindMin(RowMin, rowcount);
    max = FindMax(RowMax, rowcount);
//    printf("max value of data is %f\n", max);
    free(RowMin);
    free(RowMax);
    
    double denomintor = max-min;
    for (int i = 0; i<rowcount; i++) {
        for (int j = 0; j<colcount; j++) {
            data[i][j] = (data[i][j]-min)/denomintor;
        }
    }

}

double FindMin(double* num, int len){
    
    double min = *num;
    for (int i=0; i<len; i++) {
        if (min>*(num+i)) {
            min = *(num+i);
        }
    }
    return min;
}

double FindMax(double* num, int len){
    double max = *num;
    for (int i=0; i<len; i++) {
        if (max<*(num+i)) {
            max = *(num+i);
        }
    }
    return max;
}
