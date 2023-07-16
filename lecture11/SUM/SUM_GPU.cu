#include "SUM_GPU.cuh"
int main(int argc, char* argv[]){
    cudaStream_t* stream;
    stream=(cudaStream_t*)malloc(sizeof(cudaStream_t));
    memset( stream, 0, sizeof(cudaStream_t) );
    cudaStreamCreate(&stream[0]);
    LARGE_INTEGER CPU_Start,CPU_End;
    long long until = 2000000000;
    long long CPU_RES;
    
    //Calcu CPU cumsumption
    QueryPerformanceCounter(&CPU_Start);
    CPU_RES = CPU_SUM(until);
    QueryPerformanceCounter(&CPU_End);
    printf("Time used in CPU: %lf, result:%lld\n",cuTime(&CPU_Start,&CPU_End),CPU_RES);
    
    //Calcu sum in GPU;
    LARGE_INTEGER CudaStart,CudaEnd; 
    long long* partial;
    long long* device_ptr;
    int TotalBlock = 1024;
    int Thread = 1024;
    cudaHostAlloc(&partial,sizeof(long long)*TotalBlock,cudaHostAllocMapped);
    cudaHostGetDevicePointer(&device_ptr,partial,0);
    long long* HostRes;
    long long* Resptr;
    cudaHostAlloc(&HostRes,sizeof(long long) * 1 , cudaHostAllocMapped);
    cudaHostGetDevicePointer(&Resptr , HostRes , 0);

    printf("NBlock\tNThread\t\tResult\t\t  Time\n");
    int shared_memory = sizeof(long long) * Thread;
    QueryPerformanceCounter(&CudaStart);

    GPU_SUM <<<TotalBlock , Thread , shared_memory, *stream >>> (until,device_ptr);
    GPU_Reduction <<<1 , TotalBlock, 0 , *stream >>> (device_ptr , Resptr);
    cudaDeviceSynchronize();
    QueryPerformanceCounter(&CudaEnd);
    
    printf("%d\t%d\t%lld\t%lf\n",TotalBlock,Thread,HostRes[0],cuTime(&CudaStart,&CudaEnd));

    fflush(stdout);
    cudaFree(partial);
    cudaFree(HostRes);
    cudaFree(device_ptr);
    cudaFree(Resptr);
    return 0;
}

