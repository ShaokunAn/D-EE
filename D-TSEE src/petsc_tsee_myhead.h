//
//  petsc_tsee_myhead.h
//  petsc_tsee
//
//  Created by 安少坤 on 2019/12/9.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#ifndef petsc_tsee_myhead_h
#define petsc_tsee_myhead_h
PetscErrorCode ReadCsvData(double** data, FILE* fp);
void WriteCsvData(char *filename, double* a,int m,int n);
PetscErrorCode GetCsvRowCount(FILE *fp, PetscInt *count);
PetscErrorCode GetCsvColCount(FILE* fp, PetscInt *count);
PetscErrorCode tsee_error(Mat XX, Mat Wp, Mat Wn, PetscScalar lambda, Mat ker, PetscInt cellcount, PetscInt d, PetscScalar* e);
PetscErrorCode SolveDir(KSP ksp, Mat G, Mat P);
PetscErrorCode tseels(Mat XXold, Mat XX, Mat Wp, Mat Wn, Mat ker, PetscScalar lambda, Mat P, PetscScalar ff, Mat G, PetscScalar *alpha, PetscInt cellcount, PetscInt d);
#endif /* petsc_tsee_myhead_h */

