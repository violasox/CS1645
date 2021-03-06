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
    int n = 128;
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
    
    // Binary tree summing of all processes' integrals
    double numLevels = log2((double) numP);
    double data;
    // For each level of the tree, send and receive from appropriate processes
    for (int level = 0; level < numLevels; level++) {
        int num = pow(2, level);
        // Node receives data if it meets this condition
        if (myRank % (int) pow(2, level+1) == 0) {
            int sendingNode = myRank + num;
            if (sendingNode < numP) {
                // printf("Process %d receiving data from process %d\n", myRank, sendingNode);
                MPI_Recv(&data, 1, MPI_DOUBLE, myRank + num, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                myIntegral += data;
            }
        }
        // Node sends data if it meets this condition
        if (myRank % (int) pow(2, level+1) == num) {
            int receivingNode = myRank - num;
            // printf("Process %d sending data to process %d\n", myRank, receivingNode);
            MPI_Send(&myIntegral, 1, MPI_DOUBLE, receivingNode, 0, MPI_COMM_WORLD);
        }
    }

    double myElapsed = MPI_Wtime() - myStart;
    double totalElapsed;
    MPI_Reduce(&myElapsed, &totalElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Root reports the final integral (sum of all nodes)
    if (myRank == 0)
        printf("The integral is %f, calculated in %f seconds\n", myIntegral, totalElapsed);
    

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
