# include <math.h>
#include <time.h>
#include <stdio.h>
#include<iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<malloc.h>
#include <stdlib.h>
void sum_cpu(const int m)//纯cpu计算整数求和1+2+...+5000000000,为了体现gpu并行加速的效果，改大了数据
{
	double s = 0.0;
	for (int i = 0; i <m+1; i++)
	{
		s = s + i;
	}
	printf("s=%20.16e \n", s);
}
__global__ void sum_gpu(const int n, double* sum_local_gpu)
{
	double s = 0.0;
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	for (int j = i; j < n+1; j=j+gridDim.x*blockDim.x)
	{
		s = s + j;
	}
	sum_local_gpu[i] = s;
}
int main()
{
	const long end = 5000000000, p = 100, mp = 100;   // p 线程块数目,mp 每块线程数,ebd为求和末尾数字
	int m = p * mp;//m为总线程数
	clock_t time1, time2;
	time1 = clock();
	sum_cpu(end);
	time2 = clock();
	printf("Time for CPU run  is: %f seconds \n", (double)(time2 - time1) / CLOCKS_PER_SEC);//这里以上是cpu计算部分

	time1 = clock();
	double* sum_local_gpu;
	int MM = m * sizeof(double);//需要开辟的内存空间大小
	double* sum_local_cpu = (double*)malloc(MM);//cpu开辟内存
	cudaMalloc((void**)&sum_local_gpu, MM);//gpu开辟内存

	sum_gpu << <p, mp >> > (end, sum_local_gpu);
	cudaDeviceSynchronize();
	cudaMemcpy(sum_local_cpu, sum_local_gpu, MM, cudaMemcpyDeviceToHost);
	double sum_s = 0.0;
	for (int k = 0; k < m; ++k)//对每一个线程的求和结果再求和
	{
		sum_s += sum_local_cpu[k];
	}
	printf("sum gpu=%20.16e \n", sum_s);
	time2 = clock();
	printf("Time for GPU run  is: %f seconds \n", (double)(time2 - time1) / CLOCKS_PER_SEC);
	free(sum_local_cpu);
	cudaFree(sum_local_gpu);
	return 0;
}
