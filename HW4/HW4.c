#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "HW2.h"

// Threads per block
#define M 512
// Points in discrete integral
#define N 128
#define START_VAL 0.0
#define END_VAL 1.0

int main() {
    int numP;
    double deltaX = (END_VAL - START_VAL) / (double) N;
    int numBlocks = (N - 1)/M + 1;
    double* resultVec = double[numBlocks];
    double* deltaXDev;
    double* resultDev;
    cudaMalloc((void**)&deltaXDev, sizeof(double));
    cudaMalloc((void**)&resultDev, sizeof(double)*numBlocks);
    cudaMemcpy(deltaXDev, &deltaX, sizeof(double), cudaMemcpyHostToDevice);
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCeate(&stop);
    cudaEventRector(start);

    deviceIntegrate<<<numBlocks, M>>>(deltaXDev, resultDev);
    cudaMemcpy(&resultVec, resultDev, sizeof(double)*numBlocks, cudaMemcpyDeviceToHost);
    double integral = 0;
    for (int i = 0; i < numBlocks; i++)
        integral += resultVec[i];

    cudaEventSynchronize(stop);
    float elapsedMili = 0;
    cudaEventElapsedTime(&elapsedMili, start, stop);

    printf("The integral is %f, calculated in %f seconds\n", integral, elapsedMili / 1000.0);
    

    return(0);
}

__global__ void deviceIntegrate(double* deltaX, double* result) {
    __shared__ double resultVector[M];
    result[blockIdx.x] = 0;
    int globalIdx = threadIdx.x + blockIdx.x*M;
    if (globalIdx < N) {
        int localIdx = threadIdx.x;
        double myVal = START_VAL + globalIdx*(*deltaX);
        resultVector[localIdx] = calcStep(myVal, *deltaX);
    }
    atomicAdd(result[blockIdx.x], resultVector[localIdx]);
}

// Calulate one discrete step of the integral
double calcStep(double x, double deltaX) {
    double numerator = function(x) + function(x + deltaX);
    // printf("Process %d: Integral from %f to %f is %f\n", myRank, x, x+deltaX, (numerator*deltaX/2));
    return((numerator / 2) * deltaX);
}


// The function to calculate the integral of
// Can be changed
double function(double x) {
    return(4 / (1 + (x*x)));
}
