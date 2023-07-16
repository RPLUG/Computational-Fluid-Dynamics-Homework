#include <Windows.h>
#include <math.h>   
#include <cuda_runtime.h>
#include <profileapi.h>
#include <winnt.h>
#include <stdio.h>

//kernal funcation used in GPU
__global__ void 
GPU_SUM(const long long until,long long* device_ptr){
    extern __shared__ long long part[];
    long long sum = 0;
    const int curThread = threadIdx.x;
    for( int index = blockIdx.x * blockDim.x + curThread;
            index <= until;
            index += blockDim.x * gridDim.x)
    {
        sum += (index);
    }
    part[curThread] = sum;
    __syncthreads();
    for(int activeThread = blockDim.x >> 1;
            activeThread;
            activeThread >>= 1)
    {
        if(curThread<activeThread){
            part[curThread] += part[curThread+activeThread];
        }
        __syncthreads();
    }
    if(curThread == 0){
        device_ptr[blockIdx.x] = part[0];
    }
}
__global__ void
GPU_Reduction(long long* source, long long* ptr){
    for(int activeThread = blockDim.x >> 1;
            activeThread;
            activeThread >>= 1)
    {
        if( threadIdx.x < activeThread){
            source[threadIdx.x]+=source[threadIdx.x+activeThread];
        }
        __syncthreads();
    }
    if(threadIdx.x == 0){
        *ptr=source[0];
    }
}

//used to calcu Sum in cpu
long long CPU_SUM(long long until){
    long long res=0;
    for(int i=1;i<=until;i++){
        res+=i;
    }
    return res;
}

//Used to calculate time in host using high performance timerecord
inline double cuTime(LARGE_INTEGER* Start,LARGE_INTEGER* End){
    double diff=(double)(End->QuadPart-Start->QuadPart);
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return diff/(double) freq.QuadPart;
}
