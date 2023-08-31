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
#include <petsc.h>


PetscInt **int_get_mem(int rowcount, int colcount);
double **double_get_mem(int rowcount, int colcount);
PetscErrorCode sqdist(double* x, double* y, int len, double *dist);
PetscErrorCode merge(double* arr, PetscInt *order, int l, int m, int r);
PetscErrorCode mergeSort(double* arr, PetscInt* order, int l, int r);
PetscErrorCode nndist(double** data, double** D2, PetscInt** nn, PetscInt cellcount, PetscInt rowcount, PetscInt featcount, PetscInt k, PetscInt istart, PetscInt iend);
PetscErrorCode find_num(double* arr, double eps, PetscInt totalnum, int *num);
PetscErrorCode find(int *ind,  double* arr, double eps, PetscInt totalnum);
PetscErrorCode eabounds(double** Bbounds, double** D2, double logK, PetscInt k, PetscInt rowcount);
PetscErrorCode eabeta(double *dis, double b0, double logK, double *B, double *Wp, PetscInt k, double *out);

PetscErrorCode ea(double** data, PetscInt K, PetscInt cellcount, PetscInt rowcount, PetscInt featcount, Mat Wp){
    // compute knn distance matrix of data in each processor,
    // rowcount is the number of cells stored in each processor, not all cells' count
    PetscErrorCode ierr;
    PetscInt j, l;
    PetscInt k=cellcount-1;
    double** D2;
    double ** Bbounds;
    PetscInt** nn;
    D2 = double_get_mem(rowcount, k);
    nn = int_get_mem(rowcount, k);
    Bbounds = double_get_mem(rowcount, 2);
    double * tempWp;
    tempWp = (double*) malloc(sizeof(double)*k);
    PetscInt istart, iend, nlocal;
    ierr = MatGetOwnershipRange(Wp, &istart, &iend); CHKERRQ(ierr);
    ierr = nndist(data, D2, nn, cellcount, rowcount, featcount, k, istart, iend); CHKERRQ(ierr);
    
    // compute the upper and lower bounds for beta in each processor
    ierr = eabounds(Bbounds, D2, log(K), k, rowcount); CHKERRQ(ierr);


    //    compute beta in the order of bcal_order,
    //    and bcal_order is obtained by the order of the distance to its k-th nearest point for each point
    double *b=NULL, *Kdis=NULL;
    PetscInt *bcal_order=NULL;
    b = (double*)malloc(rowcount*sizeof(double));
    Kdis = (double*)malloc(rowcount*sizeof(double));
    ierr = PetscMalloc1(rowcount, &bcal_order); CHKERRQ(ierr);


    for (j = 0; j<rowcount; j++) {
        Kdis[j] = D2[j][K-1];
        bcal_order[j] = j;
    } 
    
    // the computation order is determined by the K-th distance of each processor,
    // not by all the K-th distance of all cells
    ierr = mergeSort(Kdis, bcal_order, 0, rowcount-1); CHKERRQ(ierr);
    j = bcal_order[0];
    double b0 = Bbounds[j][0]+Bbounds[j][1];
    b0 = b0/2;


    ierr = MatGetLocalSize(Wp, &nlocal, NULL); CHKERRQ(ierr);
    ierr = MatZeroEntries(Wp); CHKERRQ(ierr);
    PetscInt tt_order = j+istart; 

    for(l=0;l<nlocal; l++){

        ierr = eabeta(D2[j], b0, log(K), Bbounds[j], tempWp, k, b+j); CHKERRQ(ierr);

        ierr = MatSetValues(Wp, 1 , &tt_order, k, nn[j], tempWp, INSERT_VALUES); CHKERRQ(ierr);

        b0 = b[j];
        j = bcal_order[l+1];
        tt_order = j+istart;
    }
    ierr = MatAssemblyBegin(Wp, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Wp, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    free(D2);
    free(nn);
    free(tempWp);
    free(Bbounds);
    free(b);
    free(Kdis);
    ierr = PetscFree(bcal_order); CHKERRQ(ierr);
    return ierr;
}


PetscErrorCode nndist(double** data, double** D2, PetscInt** nn, PetscInt cellcount, PetscInt rowcount, PetscInt featcount, PetscInt k, PetscInt istart, PetscInt iend){
//    compute knn distance matrix and its corresponding index matrix
    PetscErrorCode ierr;
    double *tempdist=NULL;
    PetscInt *tempnnorder=NULL;
    tempdist = (double*)malloc(cellcount*sizeof(double));
    tempnnorder = (PetscInt*)malloc(cellcount*sizeof(PetscInt));

    int t1, i, j;
    for (i=0; i<rowcount; i++) {
        for (j=0; j<cellcount; j++) {
            ierr = sqdist(data[i+istart], data[j], featcount, tempdist+j); CHKERRQ(ierr);
            tempnnorder[j] = j;
        }
        
        ierr = mergeSort(tempdist, tempnnorder, 0, cellcount-1); CHKERRQ(ierr);
        
        for (t1 = 0; t1<k; t1++) {
            D2[i][t1] = tempdist[t1+1];
            nn[i][t1] = tempnnorder[t1+1];
        }
    }

    free(tempdist);
    free(tempnnorder);
    return ierr;
    
}


PetscErrorCode sqdist(double* x, double* y, int len, double *dist){

    PetscErrorCode ierr=0;
    double out=0;
    int i;
    double xi, yi;
    for (i=0; i<len; i++) {
        xi = *(x+i);
        yi = *(y+i);
        out = out + pow(xi-yi, 2);
    }
    *dist = out;
    return ierr;
}

PetscErrorCode mergeSort(double* arr, PetscInt* order, int l, int r){
    
    PetscErrorCode ierr=0;
    int m;
    if (l<r) {
        m = l+(r-l)/2;
        
        ierr = mergeSort(arr, order, l, m); CHKERRQ(ierr);
        ierr = mergeSort(arr, order, m+1, r); CHKERRQ(ierr);
        
        ierr = merge(arr, order, l, m, r); CHKERRQ(ierr);
    }
    return ierr;
}

PetscErrorCode merge(double* arr, PetscInt *order, int l, int m, int r){

    PetscErrorCode ierr=0;
    int n1 = m-l+1;
    int n2 = r-m;
    int i,j,k;
    
    double *L=NULL, *R=NULL;
    L = (double*)malloc(n1*sizeof(double));
    R = (double*)malloc(n2*sizeof(double));

    PetscInt *orderL=NULL, *orderR=NULL;
    orderL = (PetscInt*) malloc(n1*sizeof(PetscInt));
    orderR = (PetscInt*) malloc(n2*sizeof(PetscInt));

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
    free(L);
    free(R);
    free(orderL);
    free(orderR);
    return ierr;
}

PetscErrorCode eabounds(double** Bbounds, double** D2, double logK, PetscInt k, PetscInt rowcount){

    PetscErrorCode ierr = 0;
    double logN = log(k);
    double logNK = logN - logK;

    double *delta2=NULL;
    delta2 = (double*)malloc(rowcount*sizeof(double));

    int i, j, ind_num, flag=1;
    for (i=0; i<rowcount; i++) {
        delta2[i] = D2[i][1] - D2[i][0];
    }
    double eps = 2.2*pow(10, -16);
    ierr = find_num(delta2, eps, rowcount, &ind_num); CHKERRQ(ierr);
    i=2;
    if (ind_num>0) {
        int *ind=NULL;
        ind = (int*)malloc(k*sizeof(int));
        memset(ind, -1, sizeof(int)*k);
        
        while (ind_num>0) {
            ierr = find(ind, delta2, eps, rowcount); CHKERRQ(ierr);
            if (i>exp(logK) && flag) {
                for (j = 0; j<ind_num; j++){
                    D2[ind[j]][0] = D2[ind[j]][0]*0.99;
                }
                flag = 0;
            }
            for (j = 0; j<ind_num; j++) {
                delta2[ind[j]] = D2[ind[j]][i] - D2[ind[j]][0];
            }
            ierr = find_num(delta2, eps, rowcount, &ind_num); CHKERRQ(ierr);
            i++;
        }
        free(ind);
    }
    
    double *deltaN=NULL;
    deltaN = (double*)malloc(rowcount*sizeof(double)); 
    for (j = 0; j<rowcount; j++)
        deltaN[j] = D2[j][k-1]-D2[j][0];
    
    //    compute p1 = p1(k,logK)
    double p1;
    if(logK>log(sqrt(2*k)))
        p1 = 0.75;
    else{
        p1 = 0.25;
        double e, g;
        for (i = 1; i<100; i++) {
            e = -p1*log(p1/k)-logK;
            g = -log(p1/k)-1;
            p1 = p1-e/g;
        }
        p1 = 1-p1/2;
    }
    
    //    compute the upper and lower limit of log-beta
    double bU;
    double bL1;
    double bL2;
    double bL;
    for (j = 0; j<rowcount; j++) {
        bU = (2*log(p1*(k-1)/(1-p1)))/delta2[j];
        bL1 = (2*logNK/(1-1/k))/deltaN[j];
        bL2 = (2*sqrt(logNK))/sqrt(pow(D2[j][k-1],2)-pow(D2[j][0], 2));
        bL = bL1>bL2 ? bL1 : bL2;
        Bbounds[j][0] = log(bL);
        Bbounds[j][1] = log(bU);
    }

    free(delta2);
    free(deltaN);

    return ierr;

}


PetscErrorCode find_num(double* arr, double eps, PetscInt totalnum, int *num){
    //    only used to find the number of index of which the value is LESS to eps
    PetscErrorCode ierr=0;
    *num=0;
    int i;
    for (i=0; i<totalnum; i++) {
        if(arr[i]<eps)
            (*num)++;
    }
    return ierr;
}

PetscErrorCode find(int *ind,  double* arr, double eps, PetscInt totalnum){

    PetscErrorCode ierr=0;
    int i, k=0;
    
    for (i = 0; i<totalnum; i++) {
        if (arr[i]<eps) {
            ind[k] = i;
            k++;
        }
    }
    return ierr;
}

PetscErrorCode eabeta(double *dis, double b0, double logK, double *B, double *Wp, PetscInt k, double *out){
    
    PetscErrorCode ierr=0;
    double b = *out;
    int maxit=20, i=1, pbm;
    double tol = pow(10, -10);
    
    if(b0<B[0] || b0>B[1])
        b = (B[0]+B[1])/2;
    else
        b = b0;
    
    double dispold,dispnew, bE, m0, m1, m2, m12, e, eg2, g, esqd, p;
    double *ed2;
    ed2 = (double*) malloc(k*sizeof(double));
    
    int j;
    while (1) {
        dispold = B[1]-B[0];
        bE = exp(b);
        pbm = 0;
        m0 = 0;
        m1 = 0;
        m2 = 0;
        
        //        compute the function value:m0, m1v(m1 vevtor), m1
        
        for (j = 0; j<k; j++) {
            ed2[j] = exp(-dis[j]*bE);
            m0 = m0+ed2[j];
        }
        
        if (m0<DBL_MIN) {
            e = -logK;
            pbm = 1;
        }
        else{
            for (j=0; j<k; j++) {
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
            for (j=0; j<k; j++) {
                m2 = m2+dis[j]*ed2[j]*dis[j]/m0;
            }
            m12 = pow(m1, 2)-m2;
            g = eg2*m12;// based on the formulation of gradient, g = bE*m12, not g = bE*bE*m12, the additional bE can be regarded as step size during iteration
            
            if(g==0)
                pbm = 1;
        }
        
        if (pbm) {
            esqd = 0;
            for (j = 0; j<k; j++) {
                esqd = esqd + exp(-dis[j]* exp(B[0])) + exp(-dis[j]*exp(B[1]));
            }
            if(esqd<2*sqrt(DBL_MIN))
                break;
            b = (B[0]+B[1])/2;
            i=1;
            continue;
        }
        
        // Newton step ok, update bounds
        p = -e/g;
        b = b+p;
        
        if (b<B[0] || b>B[1]) {
            b = (B[0] + B[1])/2;
            i = 0;
        }
        i++;
    }
    
    for (j = 0; j<k; j++)
        Wp[j] = ed2[j]/m0;
    
    free(ed2);
    *out = b;
    return ierr;
}

