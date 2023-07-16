#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include<time.h>
using namespace std;

int main(int argc, char* argv[])
{
	double t1, t2;
	int N = 2000;//考虑计算机性能，缩减问题规模为2000*2000
	int myid, P;
	double A[2000][2000];
	double B[2000][2000];
	double C[2000][2000];
	double x; double y; double s = 0;	
	for (int i = 0; i < N; i++)//A、B矩阵初始化
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
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//获取进程编号
	MPI_Comm_size(MPI_COMM_WORLD, &P);//获取进程数
	MPI_Status status;
	int NP = N / P;
	double* A1 = new double [NP*N];//对每一个进程分割出A1(NP,N),B1(N,NP)和C1(NP,NP)，由于MPI通信的限制，必须使用一维数组传递信息
	for (int i = 0; i < NP; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A1[i * N + j] = A[i + myid * NP][j];//每一个进程的A1赋初值
		}
	}
	double* B1 = new double [N*NP];
	double* Btemp = new double [N*NP];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < NP; j++)
		{
			B1[i*NP+j] = B[i][j+ myid * NP];//每一个进程的B1赋初值
			//Btemp[i][j] = B[i][j+ myid * NP];
		}
	}
	//double* C1 = new double [NP];
	double** Ctemp = new double* [NP];//开辟Ctemp的存储空间，赋值为0，这里把组装myid块矩阵包含在了循环里
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
	//cout <<"test:B1第一项值" << B1[0] <<"进程id为"<<myid << endl;
	double temp; double tempmyid=0;
	double res = 0;
	for (int step = 0; step < P; step++)
	{	double tempres = 0;
		int id_send =(myid + step>P-1 ? myid + step - P: myid + step);
		int id_recv =(myid - step<0 ? myid - step + P: myid - step);
		MPI_Sendrecv(B1, N * NP, MPI_DOUBLE, id_send, 99,  Btemp, N * NP, MPI_DOUBLE, id_recv, 99,MPI_COMM_WORLD,&status);//B1传入工作空间Btemp
		for (int i = 0; i < NP; i++)
		{
			for (int j = 0; j < NP; j++)
			{	
				double sum = 0;//每次求Ctemp前对sum置0
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
						C[m][n] = Ctemp[i][j];//给C赋值，注意不能在进程里对C直接求和，因为对于每个进程C都是不完整的
					}					
				}*/
				tempres =tempres+ Ctemp[i][j] * Ctemp[i][j];//临时储存矩阵C的一块的各项平方和
			}
		}
		tempmyid = tempmyid + tempres;//储存一个进程计算的矩阵C的平方和，即C1各项平方和
	}
	//cout <<"test:tempres值" << tempmyid <<"进程id为"<<myid << endl;
	delete[] Ctemp;
	delete[] Btemp;
	delete[] B1;
	delete[] A1;
	t2 = MPI_Wtime();
	MPI_Reduce(&tempmyid, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//用规约命令把各进程结果结合起来
	res = res *1.0/ N / N;
	if (myid == 0)
	{
		cout << P << "个进程并行计算耗时" << t2 - t1 << "s" << endl;
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
		cout << "串行结果为S=" << lineres;
	}
	MPI_Finalize();
	
	return 0;
}