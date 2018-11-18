#include <math.h>
#include <pthread.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "HW3.h"

int numThreads;
double startVal;
double endVal;
int n;
double integral = 0;
pthread_mutex_t mutex;

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Arguments: number of threads, starting value, ending value, number of points\n");
        return 1;
    }
    sscanf(argv[1], "%d", &numThreads);
    sscanf(argv[2], "%lf", &startVal);
    sscanf(argv[3], "%lf", &endVal);
    sscanf(argv[4], "%d", &n);

    clock_t start = clock();

    // Check that n is evenly divisible by p
    if (n % numThreads != 0) {
        printf("Please make sure n is evenly divisible by p\n");
        return 1;
    }
    pthread_mutex_init(&mutex, NULL);
    pthread_t* threadHandles;
    threadHandles = malloc(numThreads*sizeof(pthread_t));

    for (long i = 0; i < numThreads; i++) 
        pthread_create(&threadHandles[i], NULL, work, (void*) i);
    
    for (int i = 0; i < numThreads; i++)
        pthread_join(threadHandles[i], NULL);
        
    pthread_mutex_destroy(&mutex);

    double timeElapsed = (double) (clock() - start) / CLOCKS_PER_SEC;
    // Root reports the final integral (sum of all nodes)
    printf("The integral is %f, calculated in %f seconds\n", integral, timeElapsed);
    return(0);
}

void* work(void* input) {
    long myRank = (long) input;
    double deltaX = (endVal - startVal) / n;
    // Points in the integral that each process handles
    int pointsPerNode = n / numThreads;
    // Starting point for this process
    double curPoint = startVal + myRank*pointsPerNode*deltaX;
    double myIntegral = 0;
    // Calculate the integral for this process' section
    for (int i = 0; i < pointsPerNode; i++) {
        myIntegral += calcStep(curPoint, deltaX, myRank);
        curPoint += deltaX;
    }
    pthread_mutex_lock(&mutex);
    integral += myIntegral;
    pthread_mutex_unlock(&mutex);
}

// Calulate one discrete step of the integral
double calcStep(double x, double deltaX, int myRank) {
    double numerator = function(x) + function(x + deltaX);
    // printf("Process %d: Integral from %f to %f is %f\n", myRank, x, x+deltaX, (numerator*deltaX/2));
    return((numerator / 2) * deltaX);
}


// The function to calculate the integral of
// Can be changed
double function(double x) {
    return(4 / (1 + (x*x)));
}
