//
//  EE_main.c
//
//  Created by 安少坤 on 2019/3/25.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "myhead.h"
#include <math.h>


int main(int argc, const char * argv[]) {
    
    char filename[100];
    strcpy(filename, argv[1]);

    FILE* fp;
    fp = fopen(filename, "r");
    if (fp==NULL)
    {
        printf("Error: Fail to open CSV file!\n");
        fclose(fp);
        exit(-1);
    }
    
    printf("Successfully open CSV file!\n");

    int colcount, rowcount;
    colcount = GetCsvColCount(fp);
    rowcount = GetCsvRowCount(fp);

    printf("feature count = %d\ncell count = %d\n",colcount, rowcount);
    
    //        allocate dynamic memory for data
    double** data;
    data = double_get_mem(rowcount, colcount);

    ReadCsvData(data, fp);
    
    NormalizeData(data, rowcount, colcount);
    
    //    compute the entropic affinity for each point
    double ** Wp;
    double ** Wn;
    Wp = double_get_mem(rowcount, rowcount);
    Wn = double_get_mem(rowcount, rowcount);

    
//    set the entropy of distribution for each cell to be logK, and K defaults 20
    ea(data, 20, rowcount, colcount, Wp);
    
    for (int i = 0; i<rowcount; i++) {
        for (int j=0; j<rowcount; j++) {
            Wn[i][j] = sqdist(data[i], data[j], colcount);
        }
    }
    
//    make sure both matrices are normalized, symmetric and have zeros on the diagonal
    double Wp_sum = 0, Wn_sum = 0;
    for (int i = 0; i<rowcount; i++) {
        for (int j = 0; j<i; j++) {
            Wn[i][j] = (Wn[i][j]+Wn[j][i])/2;
            Wn[j][i] = Wn[i][j];
            Wn_sum = Wn_sum+Wn[i][j];
            Wp[i][j] = (Wp[i][j]+Wp[j][i])/2;
            Wp[j][i] = Wp[i][j];
            Wp_sum = Wp_sum+Wp[i][j];
        }
        Wn[i][i] = 0;
        Wp[i][i] = 0;
    }
    Wn_sum = 2*Wn_sum;
    Wp_sum = 2*Wp_sum;
 
    for (int i=0; i<rowcount; i++) {
        for (int j=0; j<i; j++) {
            Wn[i][j] = Wn[i][j]/Wn_sum;
            Wn[j][i] = Wn[i][j];
            Wp[i][j] = Wp[i][j]/Wp_sum;
            Wp[j][i] = Wp[i][j];
        }
    }

    
    int d = 2; // the required dimension of embeddings we want
    int maxit = 100;
    double tol = pow(10, -3);
    
    double lambda = 10; // the weight to balance attractive and repulsive force terms
    
    double** XX;
    XX = double_get_mem(rowcount, 2);

    char X0Filename[100];
    strcpy(X0Filename, argv[2]);
    FILE* X0_fp;
    X0_fp = fopen(X0Filename, "r");
    if (X0_fp==NULL) {
        printf("Error: Fail to open X0 csv file!\n");
        fclose(X0_fp);
        exit(-1);
    }
    
    ReadCsvData(XX, X0_fp);
    
    
    ee(Wp, Wn, d, lambda, XX, rowcount, maxit, tol);
    
    double *x;
    x = (double*)malloc(rowcount*d*sizeof(double));
    int i, j;
    for(i=0; i<rowcount; i++){
        for(j=0; j<d; j++){
            x[i*d+j] = XX[i][j];
        }
    }
    char outname[100];
    strcpy(outname, argv[3]);

    WriteCsvData(outname, x, rowcount, d);
    
    free(x);
    free(XX);
    free(Wp);
    free(Wn);
    free(data);
    
}



