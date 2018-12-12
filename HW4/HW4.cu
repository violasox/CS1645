#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "HW4.h"

// Threads per block
#define THREADS_PER_BLOCK 32
// Points in discrete integral
#define N 1000000000
#define NUM_THREADS 1024
// Integral begins at
#define START_VAL 0.0
// Integral ends at
#define END_VAL 1.0

int main() {
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    int numBlocks = (NUM_THREADS - 1)/THREADS_PER_BLOCK + 1;
    printf("Calculating integral with %d points and %d GPU blocks\n", N, numBlocks);
    double* resultVec = (double*) malloc(sizeof(double)*numBlocks);
    double* resultDev;
    cudaMalloc((void**)&resultDev, sizeof(double)*numBlocks);
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    
    deviceIntegrate<<<numBlocks, THREADS_PER_BLOCK>>>(resultDev);
    cudaEventRecord(stop);
    cudaMemcpy(resultVec, resultDev, sizeof(double)*numBlocks, cudaMemcpyDeviceToHost);
    double integral = 0;
    for (int i = 0; i < numBlocks; i++) {
        integral += (double) resultVec[i];
    }

    cudaEventSynchronize(stop);
    float elapsedMili = 0;
    cudaEventElapsedTime(&elapsedMili, start, stop);

    printf("The integral is %f, calculated in %f seconds\n", integral, elapsedMili / 1000.0);
    

    return(0);
}

__global__ void deviceIntegrate(double* result) {
    __shared__ double resultVector[THREADS_PER_BLOCK];
    double deltaX = (END_VAL - START_VAL) / (double) N;
    int pointsPerThread = (N - 1)/NUM_THREADS + 1;
    result[blockIdx.x] = 0;
    int globalIdx = threadIdx.x + blockIdx.x*THREADS_PER_BLOCK;
    int localIdx = threadIdx.x;
    resultVector[localIdx] = 0;
    for (int i = globalIdx*pointsPerThread; i < (globalIdx+1)*pointsPerThread && i < N; i++) {
        double myVal = START_VAL + i*deltaX;
        resultVector[localIdx] += calcStep(myVal, deltaX);
    }
    // result[0] = 5;
    atomicAdd(&result[blockIdx.x], resultVector[localIdx]);
}

// Calulate one discrete step of the integral
__device__ double calcStep(double x, double deltaX) {
    double numerator = function(x) + function(x + deltaX);
    // printf("Process %d: Integral from %f to %f is %f\n", myRank, x, x+deltaX, (numerator*deltaX/2));
    return((numerator / 2) * deltaX);
}


// The function to calculate the integral of
// Can be changed
__device__ double function(double x) {
    return(4 / (1 + (x*x)));
}
