#include "mpi.h"
#include <stdio.h>
#include <math.h>
using namespace std;
int sum(int& end) 
{
	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//��ȡ���̱��
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//��ȡ������
	int temp = 0;//�洢һ�������ڲ�ѭ������ʱ����
	int res = 0;//�洢n�����̵ļ�����
	for (int i = myid; i <= end; i = i + numprocs) 
	{
		temp = temp + i;
	}//�ֳ�n�����̣�ÿ�����̴ӱ�ſ�ʼ�ӣ�����n=4����0+4+8...��1+5+9......
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