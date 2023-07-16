#include <Windows.h>
#include <math.h>   
#include <cuda_runtime.h>
#include <profileapi.h>
#include <winnt.h>
#include <stdio.h>

__global__ void 
CUDA_Matrix_Mult(int Total_Size,double* A,double* B,double* Out,int* Matrix_Size){
    __shared__ double storage;
    const int Blc_ID = blockIdx.x;
    const int Thr_ID = threadIdx.x;
    for(int i = 0 ; i < Matrix_Size[0] ; i ++){
        storage += A[Blc_ID*Matrix_Size[0] + i] * B[Thr_ID*Matrix_Size[1]+i];
    }
    Out[Blc_ID*Matrix_Size[0]+Thr_ID]=storage;
}

//Time Record
inline double cuTime(LARGE_INTEGER* Start,LARGE_INTEGER* End){
    double diff=(double)(End->QuadPart-Start->QuadPart);
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return diff/(double) freq.QuadPart;
}
//Used in cpu, A Stored in cols,B Stored in rows,C Stored in cols;
void Matrix_test(double* A , double* B , double* C,int* Matrix_size){
    for(int i=0 ; i<Matrix_size[0];i++){
        for(int k=0;k<Matrix_size[0];k++){
            C[i*Matrix_size[0]+k]+=A[i*Matrix_size[0]+k]*B[i*Matrix_size[0]+k];
        }
    }
}
double Matrix_SUM(double* C,int Total_Size){
    double res=0;
    for(int i=0;i<Total_Size;i++){
        res+=C[i];
    }
    return res;
}