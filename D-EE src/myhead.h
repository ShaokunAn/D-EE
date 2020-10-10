//
//  myhead.h
//  parallel_test
//
//  Created by 安少坤 on 2019/3/26.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#ifndef myhead_h
#define myhead_h

int **int_get_mem(PetscInt rowcount, PetscInt colcount);
double **double_get_mem(PetscInt rowcount, PetscInt colcount);
PetscErrorCode GetCsvRowCount(FILE *fp, PetscInt *count);
PetscErrorCode GetCsvColCount(FILE* fp, PetscInt *count);
PetscErrorCode ReadCsvData(double** data, FILE* fp);
PetscErrorCode NormalizeData(double** data, PetscInt rowcount, PetscInt colcount);
PetscErrorCode WriteCsvData(char *filename, PetscScalar *a,PetscInt m,PetscInt n);
PetscErrorCode ea(double** data, PetscInt K, PetscInt cellcount, PetscInt rowcount, PetscInt featcount, Mat Wp);
PetscErrorCode sqdist(double* x, double* y, int len, double *dist);
PetscErrorCode ee(Mat Wp, Mat Wn, PetscInt d, PetscScalar lambda, Mat XX, PetscInt allrowcount, PetscInt rowcount, PetscInt maxit, PetscScalar tol);


#endif /* myhead_h */
