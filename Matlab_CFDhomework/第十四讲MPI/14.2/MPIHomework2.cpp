#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include<time.h>
using namespace std;

int main(int argc, char* argv[])
{
	double t1, t2;
	int N = 2000;//���Ǽ�������ܣ����������ģΪ2000*2000
	int myid, P;
	double A[2000][2000];
	double B[2000][2000];
	double C[2000][2000];
	double x; double y; double s = 0;	
	for (int i = 0; i < N; i++)//A��B�����ʼ��
	{
		for (int j = 0; j < N; j++)
		{
			x = i * 1.0 / (N - 1);
			y = j * 1.0 / (N - 1);
			A[i][j] = exp(y) * sin(3.0 * x);
			B[i][j] = (x + cos(4.0 * x)) * (1.0 + y);
		}
	}
	MPI_Init(&argc, &argv);
	t1 = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//��ȡ���̱��
	MPI_Comm_size(MPI_COMM_WORLD, &P);//��ȡ������
	MPI_Status status;
	int NP = N / P;
	double* A1 = new double [NP*N];//��ÿһ�����̷ָ��A1(NP,N),B1(N,NP)��C1(NP,NP)������MPIͨ�ŵ����ƣ�����ʹ��һά���鴫����Ϣ
	for (int i = 0; i < NP; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A1[i * N + j] = A[i + myid * NP][j];//ÿһ�����̵�A1����ֵ
		}
	}
	double* B1 = new double [N*NP];
	double* Btemp = new double [N*NP];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < NP; j++)
		{
			B1[i*NP+j] = B[i][j+ myid * NP];//ÿһ�����̵�B1����ֵ
			//Btemp[i][j] = B[i][j+ myid * NP];
		}
	}
	//double* C1 = new double [NP];
	double** Ctemp = new double* [NP];//����Ctemp�Ĵ洢�ռ䣬��ֵΪ0���������װmyid������������ѭ����
	for (int i = 0; i < NP; i++)
	{
		//C1[i] = new double[NP];
		Ctemp[i] = new double[NP];
		for (int j = 0; j < NP; j++)
		{
			//C1[i][j] = 0;
			Ctemp[i][j] = 0;
		}
	}
	//cout <<"test:B1��һ��ֵ" << B1[0] <<"����idΪ"<<myid << endl;
	double temp; double tempmyid=0;
	double res = 0;
	for (int step = 0; step < P; step++)
	{	double tempres = 0;
		int id_send =(myid + step>P-1 ? myid + step - P: myid + step);
		int id_recv =(myid - step<0 ? myid - step + P: myid - step);
		MPI_Sendrecv(B1, N * NP, MPI_DOUBLE, id_send, 99,  Btemp, N * NP, MPI_DOUBLE, id_recv, 99,MPI_COMM_WORLD,&status);//B1���빤���ռ�Btemp
		for (int i = 0; i < NP; i++)
		{
			for (int j = 0; j < NP; j++)
			{	
				double sum = 0;//ÿ����Ctempǰ��sum��0
				for (int k = 0; k < N; k++)
				{
					temp = A1[i*N+k] * Btemp[k*NP+j];
					sum = sum + temp;
				}
				Ctemp[i][j] = sum;	
				
				/*for (int m = myid * NP; m < myid * (NP + 1); m++)
				{
					for (int n = step * NP; n < step * (NP + 1); n++)
					{
						C[m][n] = Ctemp[i][j];//��C��ֵ��ע�ⲻ���ڽ������Cֱ����ͣ���Ϊ����ÿ������C���ǲ�������
					}					
				}*/
				tempres =tempres+ Ctemp[i][j] * Ctemp[i][j];//��ʱ�������C��һ��ĸ���ƽ����
			}
		}
		tempmyid = tempmyid + tempres;//����һ�����̼���ľ���C��ƽ���ͣ���C1����ƽ����
	}
	//cout <<"test:tempresֵ" << tempmyid <<"����idΪ"<<myid << endl;
	delete[] Ctemp;
	delete[] Btemp;
	delete[] B1;
	delete[] A1;
	t2 = MPI_Wtime();
	MPI_Reduce(&tempmyid, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//�ù�Լ����Ѹ����̽���������
	res = res *1.0/ N / N;
	if (myid == 0)
	{
		cout << P << "�����̲��м����ʱ" << t2 - t1 << "s" << endl;
		cout << "S=" << res<<endl; 
		double lineres = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				double Ctempline = 0;
				for (int k = 0; k < N; k++)
				{
					Ctempline = Ctempline + A[i][k] * B[k][j];
				}
				lineres = Ctempline * Ctempline + lineres;
			}
		}
		lineres = lineres * 1.0 / N / N;
		cout << "���н��ΪS=" << lineres;
	}
	MPI_Finalize();
	
	return 0;
}