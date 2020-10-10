//
//  ea.c
//  parallel_test
//
//  Created by 安少坤 on 2019/3/26.
//  Copyright © 2019 Shaokun An. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

int **int_get_mem(int rowcount, int colcount);
double **double_get_mem(int rowcount, int colcount);
double sqdist(double* x, double* y, int len);
void merge(double* arr, int * order,  int l, int m, int r);
void mergeSort(double* arr, int * order, int l, int r);
void nndist(double** data, double** D2, int** nn, int cellcount, int featcount, int k);
int find_num(double*, double, int);
void find(int *ind,  double* arr, double eps, int totalnum);
void eabounds(double** Bbounds, double** D2, double logK, int k, int cellcount);
double eabeta(double * dis, double b0, double logK, double * B, double * Wp, int k);
void WriteCsvDataint_temp(char *filename, int *a,int m,int n);
void WriteCsvDatadouble_temp(char *filename, double *a,int m,int n);




void ea(double** data, int K, int cellcount, int featcount, double** Wp){
//    compute knn distance matrix
    int k = cellcount-1;
    double** D2;
    double ** Bbounds;
    int** nn;
    D2 = double_get_mem(cellcount, k);
    nn = int_get_mem(cellcount, k);
    Bbounds = double_get_mem(cellcount, 2);
    
    double * tempWp;
    tempWp = (double*) malloc(sizeof(double)*k);
    if (tempWp == NULL) {
        printf("Fail to allocate  memoey for Wp!\n");
        exit(-1);
    }
    
    nndist(data, D2, nn, cellcount, featcount, k); // compute the knn distance matrix and knn order matrix of data
    // int *nn_temp;
    // double *D2_temp;
    // nn_temp = (int*)malloc(cellcount*k*sizeof(int));
    // D2_temp = (double*)malloc(cellcount*k*sizeof(double));
    // for(int i=0; i<cellcount; i++){
    //     for(int j=0; j<k; j++){
    //         nn_temp[i*k+j] = nn[i][j];
    //         D2_temp[i*k+j] = D2[i][j];
    //     }
            
    // }
    // char nn_name[100] = "paul_petsc_nn";
    // char D2_name[100] = "paul_petsc_nnD2";
    // WriteCsvDataint_temp(nn_name, nn_temp, cellcount, k);
    // WriteCsvDatadouble_temp(D2_name, D2_temp, cellcount, k);
    // free(nn_temp);
    // free(D2_temp);
//    compute the upper and lower bounds for beta

    eabounds(Bbounds, D2, log(K), k, cellcount);
    
//    compute beta in the order of bcal_order,
//    and bcal_order is obtained by the order of the distance to its k-th nearest point for each point
    double b[cellcount];
    double Kdis[cellcount];
    int bcal_order[cellcount];
    for (int j = 0; j<cellcount; j++) {
        Kdis[j] = D2[j][K-1];
        bcal_order[j] = j;
    }
    
    mergeSort(Kdis, bcal_order, 0, cellcount-1);
    int j = bcal_order[0];
    double b0 = Bbounds[j][0]+Bbounds[j][1];
    b0 = b0/2;

    for(int l = 0; l<cellcount; l++){
        b[j] = eabeta(D2[j], b0, log(K), Bbounds[j], tempWp, k);
        // printf("%d\t %d\t %.16f\n", l, j, b0);
        memset(Wp[j], 0, cellcount);
        for(int h = 0; h<k; h++){
            Wp[j][nn[j][h]] = tempWp[h];
        }
        b0 = b[j];
        j = bcal_order[l+1];
    }
    
    free(D2);
    free(nn);
    free(tempWp);
    free(Bbounds);
}


void nndist(double** data, double** D2, int** nn, int cellcount, int featcount, int k){
//    compute knn distance matrix and its corresponding index matrix
    double tempdist[cellcount];
    int tempnnorder[cellcount];
    
    int t1;
    for (int i=0; i<cellcount; i++) {

        for (int j=0; j<cellcount; j++) {
            tempdist[j] = sqdist(data[i], data[j], featcount);
            tempnnorder[j] = j;
        }
        
        mergeSort(tempdist, tempnnorder, 0, cellcount-1);
        
        for (t1 = 0; t1<k; t1++) {
            D2[i][t1] = tempdist[t1+1];

            nn[i][t1] = tempnnorder[t1+1];
        }
    }

}


double sqdist(double* x, double* y, int len){
    double dist = 0;
    for (int i=0; i<len; i++) {
        dist = dist + pow(*(x+i)-*(y+i), 2);
    }
    return dist;
}

void mergeSort(double* arr, int* order, int l, int r){
    
    if (l<r) {
        int m = l+(r-l)/2;
        
        mergeSort(arr, order, l, m);
        mergeSort(arr, order, m+1, r);
        
        merge(arr, order, l, m, r);
    }
}

void merge(double* arr, int * order, int l, int m, int r){
    int n1 = m-l+1;
    int n2 = r-m;
    int i,j,k;
    
    double L[n1], R[n2];
    int orderL[n1], orderR[n2];
    for (i=0; i<n1; i++) {
        L[i] = arr[l+i];
        orderL[i] = order[l+i];
    }
    for (j=0; j<n2; j++) {
        R[j] = arr[j+m+1];
        orderR[j] = order[j+m+1];
    }
    
    i=0;
    j=0;
    k=l;
    while (i<n1 && j<n2) {
        if (L[i]<=R[j]) {
            arr[k] = L[i];
            order[k] = orderL[i];
            i++;
        }
        else{
            arr[k] = R[j];
            order[k] = orderR[j];
            j++;
        }
        k++;
    }
    while (i<n1) {
        arr[k] = L[i];
        order[k] = orderL[i];
        i++;
        k++;
    }
    
    while (j<n2) {
        arr[k] = R[j];
        order[k] = orderR[j];
        j++;
        k++;
    }
}


void eabounds(double** Bbounds, double** D2, double logK, int k, int cellcount){
    double logN = log(k);
    double logNK = logN - logK;
    
    double delta2[cellcount];
    for (int i=0; i<cellcount; i++) {
        delta2[i] = D2[i][1] - D2[i][0];
    }
    int ind_num;
    double eps = 2.2*pow(10, -16);
    ind_num = find_num(delta2, eps, cellcount);
    int i=2;
    int flag = 1;
    if (ind_num>0) {
        int ind[k];
        memset(ind, -1, sizeof(ind));
        
        while (ind_num>0) {
            find(ind, delta2, eps, cellcount);
            if (i>exp(logK) && flag) {
                for (int j = 0; j<ind_num; j++){
                    D2[ind[j]][0] = D2[ind[j]][0]*0.99;
                }
                flag = 0;
            }
            for (int j = 0; j<ind_num; j++) {
                delta2[ind[j]] = D2[ind[j]][i] - D2[ind[j]][0];
            }
            ind_num = find_num(delta2, eps, cellcount);
            i++;
        }
    }
    
    double deltaN[cellcount];
    for (int j = 0; j<cellcount; j++)
        deltaN[j] = D2[j][k-1]-D2[j][0];
    
//    compute p1 = p1(k,logK)
    double p1;
    if(logK>log(sqrt(2*k)))
        p1 = 0.75;
    else{
        p1 = 0.25;
        double e, g;
        for (int i = 1; i<100; i++) {
            e = -p1*log(p1/k)-logK;
            g = -log(p1/k)-1;
            p1 = p1-e/g;
//            printf("%f, %f, %f\n", e, g, p1);
        }
        p1 = 1-p1/2;
    }
    
//    compute the upper and lower limit of log-beta
    double bU;
    double bL1;
    double bL2;
    double bL;
    for (int j = 0; j<cellcount; j++) {
        bU = (2*log(p1*(k-1)/(1-p1)))/delta2[j];
        bL1 = (2*logNK/(1-1/k))/deltaN[j];
        bL2 = (2*sqrt(logNK))/sqrt(pow(D2[j][k-1],2)-pow(D2[j][0], 2));
        bL = bL1>bL2 ? bL1 : bL2;
        Bbounds[j][0] = log(bL);
        Bbounds[j][1] = log(bU);
    }

    
}

int find_num(double* arr, double eps, int totalnum ){
    //    only used to find the number of index of which the value is LESS to eps
    int num=0;
    for (int i=0; i<totalnum; i++) {
        if(arr[i]<eps)
            num++;
    }
    return num;
}

void find(int *ind,  double* arr, double eps, int totalnum){
    int k = 0;
    
    for (int i = 0; i<totalnum; i++) {
        if (arr[i]<eps) {
            ind[k] = i;
            k++;
        }
    }
}

double eabeta(double * dis, double b0, double logK, double * B, double * Wp, int k){
    double b;
    int maxit = 20;
    double tol = pow(10, -10);
    
    if(b0<B[0] || b0>B[1])
        b = (B[0]+B[1])/2;
    else
        b = b0;
    
    int i = 1;
    int pbm;
    double dispold,dispnew, bE, m0, m1, m2, m12, e, eg2, g, esqd, p;
    double ed2[k];
    
    while (1) {
        dispold = B[1]-B[0];
        bE = exp(b);
        pbm = 0;
        m0 = 0;
        m1 = 0;
        m2 = 0;
        
        //        compute the function value:m0, m1v(m1 vevtor), m1
        
        for (int j = 0; j<k; j++) {
            ed2[j] = exp(-dis[j]*bE);
            m0 = m0+ed2[j];
        }
        
        if (m0<DBL_MIN) {
            e = -logK;
            pbm = 1;
        }
        else{
            for (int j=0; j<k; j++) {
                m1 = m1+ed2[j]*dis[j]/m0;
            }
            e = bE*m1+log(m0)-logK;
        }
        
        if(fabs(e)<tol)
            break;
        
        //        very narrow bounds, no need to iterate. This can happen if K is very small
        if(B[1]-B[0]<2.22*pow(10, -16))
            break;
        
        //        update the bounds
        if (e<0 && b<=B[1])
            B[1] = b;
        else if (e>0 && b>=B[0])
            B[0] = b;
        
        dispnew = B[1] - B[0];
        if(dispold==dispnew)
            break;
        
        pbm = pbm || isnormal(e) || (e<-logK) || (e>log(k)-logK);
        
        if (!pbm){
            if(i==maxit){
                b = (B[0]+B[1])/2;
                i = 1;
                continue;
            }
            eg2 = pow(bE, 2);
            for (int j=0; j<k; j++) {
                m2 = m2+dis[j]*ed2[j]*dis[j]/m0;
            }
            m12 = pow(m1, 2)-m2;
            g = eg2*m12;// based on the formulation of gradient, g = bE*m12, not g = bE*bE*m12, the additional bE can be regarded as step size during iteration
            
            if(g==0)
                pbm = 1;
        }
        
        if (pbm) {
            esqd = 0;
            for (int j = 0; j<k; j++) {
                esqd = esqd + exp(-dis[j]* exp(B[0])) + exp(-dis[j]*exp(B[1]));
            }
            if(esqd<2*sqrt(DBL_MIN))
                break;
            b = (B[0]+B[1])/2;
            i=1;
            continue;
        }
        
        //        Newton step ok, update bounds
        p = -e/g;
        b = b+p;
        
        if (b<B[0] || b>B[1]) {
            b = (B[0] + B[1])/2;
            i = 0;
        }
        i++;
    }
    
    for (int j = 0; j<k; j++)
        Wp[j] = ed2[j]/m0;
    
    return b;
}

void WriteCsvDataint_temp(char *filename, int *a,int m,int n){
 
    // printf("\n Creating %s.csv file",filename);
    FILE *fp;
    int i,j;
    filename=strcat(filename,".csv");
    fp=fopen(filename,"w+");
    fprintf(fp, "%d", *a);
    for(i=0;i<m;i++){
        // fprintf(fp,"\n%d",i+1);
        for(j=1;j<n;j++)
            fprintf(fp,",%d",a[i*n+j]);
        if(i<m-1)
            fprintf(fp, "\n%d", a[(i+1)*n]);
    }   
    fclose(fp);
    printf("\n %s file created\n",filename);
    
}

void WriteCsvDatadouble_temp(char *filename, double *a,int m,int n){
 
    // printf("\n Creating %s.csv file",filename);
    FILE *fp;
    int i,j;
    filename=strcat(filename,".csv");
    fp=fopen(filename,"w+");
    fprintf(fp, "%.16f", *a);
    for(i=0;i<m;i++){
        // fprintf(fp,"\n%d",i+1);
        for(j=1;j<n;j++)
            fprintf(fp,",%.16f",a[i*n+j]);
        if(i<m-1)
            fprintf(fp, "\n%.16f", a[(i+1)*n]);
    }   
    fclose(fp);
    printf("\n %s file created\n",filename);
    
}
