//
//  myhead.h
//  parallel_test
//
//  Created by 安少坤 on 2019/3/26.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#ifndef myhead_h
#define myhead_h

int **int_get_mem(int rowcount, int colcount);
double **double_get_mem(int rowcount, int colcount);
int GetCsvRowCount(FILE* fp);
int GetCsvColCount(FILE* fp);
void ReadCsvData(double** data, FILE* fp);
void WriteCsvData(char *filename, double *a,int m,int n);
void NormalizeData(double** data, int rowcount, int colcount);
void ea(double** data, int K, int cellcount, int featcount, double** Wp);
double sqdist(double* x, double* y, int len);
void Cholesky(double** X, double** L, int n);
void ee(double** Wp, double** Wn, int d, double lambda, double** XX, int cellcount, int maxit, double tol);

#endif /* myhead_h */
