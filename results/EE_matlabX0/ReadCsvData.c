//
//  ReadData.c
//  parallel_test
//
//  Created by 安少坤 on 2019/3/26.
//  Copyright © 2019 Shaokun An. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_SIZE 1000000

int GetCsvRowCount(FILE *fp)
{
    //    calculate the total row count of CSV file
    int i = 0;
    char strLine[MAX_LINE_SIZE];
    fseek(fp, 0, SEEK_SET);
    while (fgets(strLine, MAX_LINE_SIZE, fp))
    {
        i++;
    }
    fseek(fp, 0, SEEK_SET);
    return i;
}

int GetCsvColCount(FILE* fp)
{
    //    calculate the total column count of CSV file
    int i = 0;
    char strLine[MAX_LINE_SIZE];
    fgets(strLine, MAX_LINE_SIZE, fp);
    fseek(fp, 0, SEEK_SET);
    char *token = strtok(strLine, ",");
    while (token) {
        i++;
        token = strtok(NULL, ",");
    }
    
    return i;
}


void ReadCsvData(double** data, FILE* fp){
    
//        Input data in CSV file
    char strLine[MAX_LINE_SIZE];
    char* str;
    int i=0;
    fseek(fp, 0, SEEK_SET);
    while (fgets(strLine, MAX_LINE_SIZE, fp)) {
        int j = 0;
        for (str = strtok(strLine, ","); str && *str; j++, str = strtok(NULL, ",")) {
            data[i][j] = atof(str);
        }
        i++;
    }
        
    fclose(fp);
}

void WriteCsvData(char *filename, double *a,int m,int n){
 
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
    
