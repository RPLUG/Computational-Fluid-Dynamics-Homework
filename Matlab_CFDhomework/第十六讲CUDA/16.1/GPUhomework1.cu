# include <math.h>
#include <time.h>
#include <stdio.h>
#include<iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<malloc.h>
#include <stdlib.h>
void sum_cpu(const int m)//��cpu�����������1+2+...+5000000000,Ϊ������gpu���м��ٵ�Ч�����Ĵ�������
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
	const long end = 5000000000, p = 100, mp = 100;   // p �߳̿���Ŀ,mp ÿ���߳���,ebdΪ���ĩβ����
	int m = p * mp;//mΪ���߳���
	clock_t time1, time2;
	time1 = clock();
	sum_cpu(end);
	time2 = clock();
	printf("Time for CPU run  is: %f seconds \n", (double)(time2 - time1) / CLOCKS_PER_SEC);//����������cpu���㲿��

	time1 = clock();
	double* sum_local_gpu;
	int MM = m * sizeof(double);//��Ҫ���ٵ��ڴ�ռ��С
	double* sum_local_cpu = (double*)malloc(MM);//cpu�����ڴ�
	cudaMalloc((void**)&sum_local_gpu, MM);//gpu�����ڴ�

	sum_gpu << <p, mp >> > (end, sum_local_gpu);
	cudaDeviceSynchronize();
	cudaMemcpy(sum_local_cpu, sum_local_gpu, MM, cudaMemcpyDeviceToHost);
	double sum_s = 0.0;
	for (int k = 0; k < m; ++k)//��ÿһ���̵߳���ͽ�������
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
