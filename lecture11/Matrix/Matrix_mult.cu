#include "Matrix_mult.cuh"
#include <random>
int main(int argc,char* argv[]){
    int N = 1024;
    int Matrix_size[2] = { N , N };
    int Total_Size = N * N * sizeof(double);
    int N_Block = 1024;
    int N_Thread = 1024;
    //----------------------HOST MALLOC--------------------------

    std::mt19937 generate{std::random_device{}()};
    std::uniform_real_distribution<double> FLY(-10000,10000);
    //A Stored in cols
    double* A=(double*)malloc(sizeof(double)*Total_Size);
    for(int i=0;i<Total_Size;i++){
        A[i]=FLY(generate);
    }
    //B Stored in rows
    double* B=(double*)malloc(sizeof(double)*Total_Size);
    for(int i=0;i<Total_Size;i++){
        B[i]=FLY(generate);
    }
    //C Stored in cols
    double* C=(double*)malloc(sizeof(double)*Total_Size);
    memset(C,0,sizeof(double)*Total_Size);

    //--------------------HOST Check------------------------------
    LARGE_INTEGER CPU_BEG , CPU_END;
    QueryPerformanceCounter( &CPU_BEG );
    Matrix_test( A , B , C , Matrix_size );
    QueryPerformanceCounter( &CPU_END );

    printf("Time Cost in the CPU are %lf, Sum of matrix are %lf\n",cuTime(&CPU_BEG,&CPU_END),Matrix_SUM(C,Total_Size));
    //--------------------CUDA MALLOC-----------------------------
    double* A_CUDA;
    double* B_CUDA;
    double* C_CUDA;
    cudaMalloc( &A_CUDA , Total_Size );
    cudaMalloc( &B_CUDA , Total_Size );
    cudaMalloc( &C_CUDA , Total_Size );
    cudaMemcpy( A_CUDA , A , Total_Size , cudaMemcpyHostToDevice );
    cudaMemcpy( B_CUDA , B , Total_Size , cudaMemcpyHostToDevice );
    cudaStream_t* stream;
    stream=(cudaStream_t*)malloc(sizeof(cudaStream_t));
    cudaStreamCreate(stream);
    LARGE_INTEGER GPU_START,GPU_END;
    //--------------------CUDA Run-----------------------------
    QueryPerformanceCounter(&GPU_START);
    CUDA_Matrix_Mult<<<N_Block,N_Thread,sizeof(double),*stream>>>(Total_Size,A_CUDA,B_CUDA,C_CUDA,Matrix_size);
    QueryPerformanceCounter(&GPU_END);
    //--------------------CUDA Check---------------------------
    cudaMemcpy(C,C_CUDA,Total_Size,cudaMemcpyDeviceToHost);
    //----------------------RESULT-----------------------------
    printf("GPU test infomation: \n");
    printf("NBlock\tNthread\t Time\t\tResult\n");
    printf("%d\t%d\t%lf\t%lf\n",N_Block,N_Thread,cuTime(&GPU_START,&GPU_END),Matrix_SUM(C,Total_Size));
    return 0;
}