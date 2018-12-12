#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "HW4.h"

// Threads per block
#define M 1024
// Points in discrete integral
#define N 50000
// Integral begins at
#define START_VAL 0.0
// Integral ends at
#define END_VAL 1.0

int main() {
    int numBlocks = (N - 1)/M + 1;
    printf("Calculating integral with %d points and %d GPU blocks\n", N, numBlocks);
    double* resultVec = (double*) malloc(sizeof(double)*numBlocks);
    double* resultDev;
    cudaMalloc((void**)&resultDev, sizeof(double)*numBlocks);
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    
    deviceIntegrate<<<numBlocks, M>>>(resultDev);
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
    __shared__ double resultVector[M];
    double deltaX = (END_VAL - START_VAL) / (double) N;
    result[blockIdx.x] = 0;
    int globalIdx = threadIdx.x + blockIdx.x*M;
    int localIdx = threadIdx.x;
    if (globalIdx < N) {
        double myVal = START_VAL + globalIdx*deltaX;
        resultVector[localIdx] = calcStep(myVal, deltaX);
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
