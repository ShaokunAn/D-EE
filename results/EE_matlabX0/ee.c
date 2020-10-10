//
//  ee.c
//
//  Created by 安少坤 on 2019/3/28.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int **int_get_mem(int rowcount, int colcount);
double **double_get_mem(int rowcount, int colcount);
double sum_prod(double* x, double* y, int start, int end);
void cholesky(double** X, double** L, int n);
double sqdist(double* x, double* y, int len);
double ee_error(double** XX, double** Wp, double** Wn, double lambda, double** ker, int cellcount, int d);
void SolveGrad(double**LTC, double** G, double** p, int n, int d);
void eels(double** XXold, double** XX, double** Wp, double** Wn, double** ker,double lambda, double** P, double *ff, double** G, double* alpha, int n, int d);
double FroNorm(double**X, int n, int d);

void ee(double** Wp, double** Wn, int d, double lambda, double** XX, int cellcount, int maxit, double tol)
{
    double Dp[cellcount];
    for(int i = 0;i<cellcount;i++){
        Dp[i] = 0;
        for(int j=0;j<cellcount;j++)
            Dp[i] = Dp[i]+Wp[i][j];
    }
    
    double** Lp4;
    Lp4 = double_get_mem(cellcount, cellcount);
       
    for (int i=0; i<cellcount; i++) {
        for (int j=0; j<i; j++) {
            Lp4[i][j] = 4*(-Wp[i][j]);
            Lp4[j][i] = Lp4[i][j];
        }
        Lp4[i][i] = 4*Dp[i];
    }
    
    //smallest nonzero element of diag(Lp4)
    double mDiagLp;
    int i=0;
    while (1 && i<cellcount) {
        if(Dp[i]>0)
            break;
        i++;
    }
    mDiagLp = Lp4[i][i];

    
    //    compute Cholesky factor of the graph Laplacian
    double OldDiagLp4[cellcount];
    for (int i = 0; i<cellcount; i++) {
        OldDiagLp4[i] = Lp4[i][i];
        Lp4[i][i] = Lp4[i][i] + pow(10, -10)*mDiagLp;
    }
    double** LTC; //L is the lower triangular cholesky factor
    LTC = (double**)malloc(cellcount*sizeof(double*));
    if (LTC==NULL) {
        printf("Fail to allocate memory for Cholesky factor!\n");
        exit(-1);
    }
    for (int i=0; i<cellcount; i++) {
        LTC[i] = (double*)malloc((i+1)*sizeof(double));
        if (LTC[i]==NULL) {
            printf("Fail to allocate memory for Cholesky factor!\n");
            exit(-1);
        }
        memset(LTC[i], 0, (i+1)*sizeof(double));
    }
    
    cholesky(Lp4, LTC, cellcount);
    
    for (int i= 0; i<cellcount; i++)
        Lp4[i][i] = OldDiagLp4[i];
    
    
//    compute low-dimension embeddings with Newton method
    double e;                               // objective value
    double** ker;//                         //ker[n,m] = exp(-|x_n -x_m|^2)
    ker = double_get_mem(cellcount, cellcount);

    e = ee_error(XX, Wp, Wn, lambda, ker, cellcount, d);
    
    double alpha = 1;
    
    int convcrit = 1;
    double** WWn;
    WWn = double_get_mem(cellcount, cellcount);
    
    
    double DDn[cellcount];
    double temp[cellcount];
    for (int i = 0; i<cellcount; i++)
        temp[i] = 1;
    
    double** G;
    double** P;
    double** L;
    double** XXold;
    double** XX_diff;
    double difffronorm, oldfronorm;
    G = double_get_mem(cellcount, d);
    P = double_get_mem(cellcount, d);
    L = double_get_mem(cellcount, cellcount);
    XXold = double_get_mem(cellcount, d);
    XX_diff = double_get_mem(cellcount, d);
    for(int i = 0; i < cellcount; i++){
        for (int j=0;j<d; j++) {
            XXold[i][j] = XX[i][j];
        }
    }
    
    
    int itecount = 1;
    while (convcrit) {
        for (int i = 0; i<cellcount; i++) {
            for (int j = 0; j<i; j++) {
                WWn[i][j] = lambda*Wn[i][j]*ker[i][j];
                WWn[j][i] = WWn[i][j];
            }
            WWn[i][i] = 0;
        }
        
        for (int i=0; i<cellcount; i++) {
            DDn[i] = sum_prod(WWn[i], temp, 0, cellcount-1);
        }
        
        for (int i = 0; i<cellcount; i++) {
            for (int j = 0; j<i; j++){
                L[i][j] = Lp4[i][j] - 4*(-WWn[i][j]);
                L[j][i] = L[i][j];
            }
            L[i][i] = Lp4[i][i] - 4*(DDn[i]-WWn[i][i]);
        }
        
        for (int i = 0; i<cellcount; i++) {
            for (int j = 0; j<d; j++) {
                G[i][j] = 0;
                for (int k = 0; k<cellcount; k++) {
                    G[i][j] += L[i][k]*XX[k][j];
                }
            }
        }
        
        SolveGrad(LTC, G, P, cellcount, d);
        
        eels(XXold, XX, Wp, Wn, ker, lambda, P, &e, G, &alpha, cellcount, d);
        e = ee_error(XX, Wp, Wn, lambda, ker, cellcount, d);
        
        for (int i = 0; i<cellcount; i++) {
            for (int j = 0; j<d; j++) {
                XX_diff[i][j] = XX[i][j] - XXold[i][j];
            }
        }
        
        difffronorm = FroNorm(XX_diff, cellcount, d);
        oldfronorm = FroNorm(XXold, cellcount, d);
        convcrit = (itecount<maxit) && (difffronorm > tol* oldfronorm);
        
        for (int i=0; i<cellcount; i++) {
            for (int j = 0; j<d; j++) {
                XXold[i][j] = XX[i][j];
            }
        }
        itecount++;
    }

    free(LTC);
    free(XX_diff);
    free(XXold);
    free(P);
    free(G);
    free(L);
    free(WWn);
    free(ker);
    free(Lp4);
}


void cholesky(double** X, double** L, int n){
    //    return the lower triangular Cholesky factor R
    L[0][0] = sqrt(X[0][0]);
    for (int i = 1; i<n; i++)
        L[i][0] = X[i][0]/L[0][0];
    
    double temp;
    for (int j = 1; j<n; j++) {
        temp = sum_prod(L[j], L[j], 0, j-1);
        L[j][j] = sqrt(X[j][j]-temp);
        
        for (int i = j+1; i<n; i++) {
            temp = sum_prod(L[i], L[j], 0, j-1);
            L[i][j] = (X[i][j] - temp)/L[j][j];
        }
        
    }
}

double sum_prod(double* x, double* y, int start, int end){
    double out = 0;
    for(int i = start;i<=end;i++)
        out += x[i]*y[i];
    return out;
}

double ee_error(double** XX, double** Wp, double** Wn, double lambda, double** ker, int cellcount, int d){
    double e = 0;
    double** sqd;
    sqd = double_get_mem(cellcount, cellcount);
    
    for (int i = 0; i<cellcount; i++) {
        for (int j=0; j<i; j++) {
            sqd[i][j] = sqdist(XX[i], XX[j], d);
            sqd[j][i] = sqd[i][j];
            ker[i][j] = exp(-sqd[i][j]);
            ker[j][i] = ker[i][j];
            e += 2*(Wp[i][j]*sqd[i][j]+lambda*Wn[i][j]*ker[i][j]);
        }
        sqd[i][i] = 0;
        ker[i][i] = 1;
        e += Wp[i][i]*sqd[i][i]+lambda*Wn[i][i]*ker[i][i];
    }
    
    free(sqd);
    return e;
}

void SolveGrad(double**LTC, double** G, double** p, int n, int d){
    //    compute -LY = G to obtain Y, then compute t(L)P = Y to obtain the P
    
    //    conpute -LY = G
    double** Y;
    Y = double_get_mem(n, d);
    
    
    double temp;
    for (int i = 0; i<d; i++) {
        Y[0][i] = -G[0][i]/LTC[0][0];
        for (int j = 1; j<n; j++) {
            temp = 0;
            for (int k = 0; k<j; k++) {
                temp += Y[k][i]*LTC[j][k];
            }
            Y[j][i] = (-G[j][i]-temp)/LTC[j][j];
        }
    }
    
    //    compute t(L)P = Y
    for (int i = 0; i<d; i++) {
        p[n-1][i] = Y[n-1][i]/LTC[n-1][n-1];
        for (int k = n-2; k>=0; k--) {
            temp = 0;
            for (int j=k+1; j<n; j++) {
                temp += p[j][i]*LTC[j][k];
            }
            p[k][i] = (Y[k][i]-temp)/LTC[k][k];
        }
    }
    
    free(Y);
}

double FroNorm(double**X, int n, int d){
    double norm = 0;
    for (int i=0; i<n; i++) {
        for (int j = 0; j<d; j++) {
            norm += pow(X[i][j], 2);
        }
    }
    norm = sqrt(norm);
    return norm;
}


void eels(double** XXold, double** XX, double** Wp, double** Wn,double** ker,double lambda, double** P, double *ff, double** G, double* alpha, int n, int d)
{
    double rho=0.8;
    double c = 0.1;
    double temp = 0;
    double e;
    double** tt;
    tt = double_get_mem(n, d);
    
    for (int i = 0; i<n ; i++) {
        for (int j = 0; j<d; j++) {
            temp += c*G[i][j]*P[i][j];
            tt[i][j] = XXold[i][j] + *alpha * P[i][j];
        }
    }
    
    
    e = ee_error(tt, Wp, Wn, lambda, ker, n, d);
    
    while (e> *ff + *alpha * temp) {
        *alpha = rho* (*alpha);
        
        for (int i=0; i<n; i++) {
            for (int j = 0; j<d; j++) {
                tt[i][j] = XX[i][j] + *alpha * P[i][j];
            }
        }
        
        e = ee_error(tt, Wp, Wn, lambda, ker, n, d);
    }
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<d; j++)
            XX[i][j] = XX[i][j]+ *alpha * P[i][j];
        
    }
    free(tt);
}
