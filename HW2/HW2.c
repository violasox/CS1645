#include <mpi.h>

int main(int argc, char** argv) {
    int numP;
    int myRank;
    double startVal = 0;
    double endVal = 1;
    int n = 10;
    double deltaX = (startVal - endVal) / (double) n;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numP);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (n % numP != 0) {
        printf("Please make sure n is evenly divisible by p\n");
        return(1);
    }

    int pointsPerNode = n / numP;
    double curPoint = startVal + myRank*pointsPerNode*deltaX;
    double myIntegral = 0;
    for (int i = 0; i < pointsPerNode; i++) {
        myIntegral += calcStep(curPoint, deltaX);
        curPoint += deltaX;
    }
    if (myRank != 0)
        MPI_Send(&myIntegral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    else {
        double data;
        for (int i = 1; i < numP; i++) {
            MPI_RECV(&data, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            myIntegral += data;
        }
        printf("The integral is %f\n", myIntegral);
    }
    return(0);
}

// Calulate one discrete step of the integral
double calcStep(double x, double deltaX) {
    double numerator = function(x) + function(x - deltaX);
    return((numerator / 2) * deltaX);
}


// The function to calculate the integral of
// Can be changed
double function(double x) {
    return(4 / (1 + (x*x)));
}
