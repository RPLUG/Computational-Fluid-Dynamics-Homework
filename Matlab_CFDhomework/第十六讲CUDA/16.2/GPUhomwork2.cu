
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<malloc.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include<cmath>

void matrix_mul_cpu(int n, float* a, float* b, float* c)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            float s = 0.0;
            for (int k = 0; k < n; ++k)
                s += a[i * n + k] * b[k * n + j];               //A(i,k)*B(k,j)
            c[i * n + j] = s;
        }
}

__global__ void matrix_mul_gpu(int n, float* a, float* b, float* c)
{
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    float s = 0.0;
    for (int k = 0; k < n; ++k)
        s += a[bid * n + k] * b[k * n + tid];               //A(i,k)*B(k,j)
    c[bid * n + tid] = s;
}
/*__global__ void transpose(int n, float* b, float* bt)
{
    bt[blockIdx.x * n + threadIdx.x] = b[threadIdx.x * n + blockIdx.x];
}
__global__ void matrix_mul_gpu2(int n, float* a, float* b, float* c)
{
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    float s = 0.0;
    for (int k = 0; k < n; ++k)
        s += a[bid * n + k] * b[tid * n + k];          //A*BT  连续访问
    c[bid * n + tid] = s;
}*/
void check_data(int n, float* c, float* c1)
{
    float s = 0.0;
    for (int i = 0; i < n * n; ++i)
        s += fabs(c[i] - c1[i]);
    printf("Total error is %f \n", s);
}
int main()
{
    const int n = 1024, M = n * n * sizeof(float);
    float* a = (float*)malloc(M);
    float* b = (float*)malloc(M);
    float* c = (float*)malloc(M);
    float* c1 = (float*)malloc(M);
    float* d_a, * d_b, * d_c, * d_bT;

    for (int i = 0; i < n * n; ++i)
    {
        a[i] = (float)(rand() % 100);
        b[i] = (float)(rand() % 100);
    }

    clock_t time1, time2, time3, time4;
    time1 = clock();
    matrix_mul_cpu(n, a, b, c);
    time2 = clock();
    printf("Time for CPU run  is: %f seconds \n", (double)(time2 - time1) / CLOCKS_PER_SEC);

    cudaMalloc(&d_a, M);
    cudaMalloc(&d_b, M);
    cudaMalloc(&d_c, M);
    cudaMemcpy(d_a, a, M, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, M, cudaMemcpyHostToDevice);

    time3 = clock();
    matrix_mul_gpu << <n, n >> > (n, d_a, d_b, d_c);  // 计算C=A*B, 每个线程计算一个点
    cudaDeviceSynchronize();
    cudaMemcpy(c1, d_c, M, cudaMemcpyDeviceToHost);
    time4 = clock();
    printf("Time for GPU run  is: %f seconds \n", (double)(time4 - time3) / CLOCKS_PER_SEC);
    check_data(n, c, c1);
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
    //cudaFree(d_bT);// by using transpose 

    /*time1 = clock();
    cudaMalloc(&d_a, M);
    cudaMalloc(&d_b, M);
    cudaMalloc(&d_c, M);
    cudaMalloc(&d_bT, M);
    cudaMemcpy(d_a, a, M, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, M, cudaMemcpyHostToDevice);

    transpose << <n, n >> > (n, d_b, d_bT);
    cudaDeviceSynchronize();
    matrix_mul_gpu2 << <n, n >> > (n, d_a, d_bT, d_c);
    cudaMemcpy(c1, d_c, M, cudaMemcpyDeviceToHost);
    time2 = clock();
    printf("Time for GPU run  is: %f seconds \n", (double)(time2 - time1) / CLOCKS_PER_SEC);
    check_data(n, c, c1);

 
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
    cudaFree(d_bT);*/ 
    free(a);
    free(b);
    free(c);
    return 0;
}




