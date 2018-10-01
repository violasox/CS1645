#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "HW2.h"

int main(int argc, char** argv) {
    int numP;
    int myRank;
    double startVal = 0;
    double endVal = 1;
    int n = 8;
    double deltaX = (endVal - startVal) / (double) n;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numP);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (n % numP != 0) {
        if (myRank == 0)
            printf("Please make sure n is evenly divisible by p\n");
        MPI_Finalize();
        return(0);
    }

    int pointsPerNode = n / numP;
    double curPoint = startVal + myRank*pointsPerNode*deltaX;
    double myIntegral = 0;
    for (int i = 0; i < pointsPerNode; i++) {
        myIntegral += calcStep(curPoint, deltaX, myRank);
        curPoint += deltaX;
    }
    
    double numLevels = log2((double) numP);
    double data;
    for (int level = 0; level < numLevels; level++) {
        int num = pow(2, level);
        if (myRank % (int) pow(2, level+1) == 0) {
            printf("Process %d receiving data from process %d\n", myRank, myRank + num);
            MPI_Recv(&data, 1, MPI_DOUBLE, myRank + num, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            myIntegral += data;
        }
        if (myRank % (int) pow(2, level+1) == num) {
            printf("Process %d sending data to process %d\n", myRank, myRank - num);
            MPI_Send(&myIntegral, 1, MPI_DOUBLE, myRank - num, 0, MPI_COMM_WORLD);
        }
    }

    if (myRank == 0)
        printf("The integral is %f\n", myIntegral);

    /*
    if (myRank != 0)
        MPI_Send(&myIntegral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    else {
        double data;
        for (int i = 1; i < numP; i++) {
            MPI_Recv(&data, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            myIntegral += data;
        }
        printf("The integral is %f\n", myIntegral);
    }
    */

    MPI_Finalize();
    return(0);
}

// Calulate one discrete step of the integral
double calcStep(double x, double deltaX, int myRank) {
    double numerator = function(x) + function(x + deltaX);
    printf("Process %d: Integral from %f to %f is %f\n", myRank, x, x+deltaX, (numerator*deltaX/2));
    return((numerator / 2) * deltaX);
}


// The function to calculate the integral of
// Can be changed
double function(double x) {
    return(4 / (1 + (x*x)));
}
