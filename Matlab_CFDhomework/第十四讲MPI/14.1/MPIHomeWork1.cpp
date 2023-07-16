#include "mpi.h"
#include <stdio.h>
#include <math.h>
using namespace std;
int sum(int& end) 
{
	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//获取进程编号
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//获取进程数
	int temp = 0;//存储一个进程内部循环的临时变量
	int res = 0;//存储n个进程的计算结果
	for (int i = myid; i <= end; i = i + numprocs) 
	{
		temp = temp + i;
	}//分成n个进程，每个进程从编号开始加，例如n=4，则0+4+8...，1+5+9......
	MPI_Reduce(&temp, &res, 1, MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD);
	return res;
}
int main(int argc, char* argv[]) {
	int numprocs;
	int myid;
	int end = 1000;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	long res = sum(end);
	if (myid == 0) {
		printf("%ld", res);
	}
	MPI_Finalize();
	return 0;
}