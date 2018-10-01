#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "HW2.h"

int main() {
    int numP;
    int myRank;
    double startVal = 0;
    double endVal = 1;
    int n = 8;
    double deltaX = (endVal - startVal) / (double) n;

    MPI_Init(NULL, NULL);
    double myStart = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &numP);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // Check that n is evenly divisible by p
    if (n % numP != 0) {
        if (myRank == 0)
            printf("Please make sure n is evenly divisible by p\n");
        MPI_Finalize();
        return(0);
    }

    // Points in the integral that each process handles
    int pointsPerNode = n / numP;
    // Starting point for this process
    double curPoint = startVal + myRank*pointsPerNode*deltaX;
    double myIntegral = 0;
    // Calculate the integral for this process' section
    for (int i = 0; i < pointsPerNode; i++) {
        myIntegral += calcStep(curPoint, deltaX, myRank);
        curPoint += deltaX;
    }
    
    double totalIntegral;
    MPI_Reduce(&myIntegral, &totalIntegral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double myElapsed = MPI_Wtime() - myStart;
    double totalElapsed;
    MPI_Reduce(&myElapsed, &totalElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Root reports the final integral (sum of all nodes)
    if (myRank == 0)
        printf("The integral is %f, calculated in %f seconds\n", totalIntegral, totalElapsed);
    

    MPI_Finalize();
    return(0);
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
