#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>
#include "solver.h"

using namespace std;

// Main function: calculates 2D Poisson with either Jacobi or Gauss-Seidel based on boolean flags and command line arguments
// First two arguments should be size of the mesh in (x y) order
// Optional 3rd argument -v triggers verbose printout mode
int main(int argc, char** argv) {
    clock_t start = clock();

    // Flag to toggle verification (true) or problem solving (false)
    bool verifyMode = true;
    
    if (argc < 5) {
        cout << "Need at least arguments: mesh size in x and y directions. Use 3rd argument -v for verbose mode.\n";
        return 1;
    }

    int totalXSize = stoi(argv[1]);
    int totalYSize = stoi(argv[2]);
    int xP = stoi(argv[3]);
    int yP = stoi(argv[4]);

    int numP;
    int myRank;
    MPI_Init(NULL, NULL);
    MPI_Comm_Size(MPI_COMM_WORLD, &numP);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if ((totalXSize - 2) % xP != 0) {
        cout << "x mesh size must be evenly divisible by the number of processors";
        return 1;
    else if ((totalYSize - 2) % yP != 0) {
        cout << "y mesh size must be evenly divisible by the number of processors";
        return 1;
    else if (numP != xP * yP) {
        cout << "Total # processors must equal x-direction processors * y-direction processors";
        return 1;
    }

    int xSize = (totalXSize - 2) / numP;
    int ySize = (totalYSize - 2) / numP;

    int xPos = myRank % xP;
    int yPos = myRank / xP;

    // Flag to trigger printing out useful info (set in command line with -v)
    bool verbose = false;
    if (argc > 5 && string(argv[5]).compare("-v") == 0)
        verbose = true;

    double xMin, xMax, yMin, yMax;
    // Set numerical mesh size
    xMin = 0;
    xMax = 2;
    yMin = 0;
    yMax = 1;
    
    double xStep = ((double) (xMax - xMin)) / ((double) (totalXSize-1));
    double yStep = ((double) (yMax - yMin)) / ((double) (totalYSize-1));
    
    const double threshold = pow(10, -12);
    
    // Allocate the mesh grid that Poisson should be calculated over
    int numElements = xSize*ySize;
    double* mesh = new double[numElements];
    setBounds(mesh, xSize, ySize, xStep, yStep, xPos, yPos);
    double* myLeftBound = new double[ySize];
    double* myRightBound = new double[ySize];
    double* myTopBound = new double[xSize];
    double* myBotBound = new double[xSize];
    double* nextLeftBound = double[ySize];
    double* nextRightBound = double[ySize];
    double* nextTopBound = double[xSize];
    double* nextBotBound = double[xSize];
    double* oldMesh = new double[numElements];
    copyMesh(mesh, oldMesh, numElements);

    if (verbose) {
        cout << "Starting mesh:\n";
        printMesh(mesh, xSize, ySize);
    }

    // Has the change after each iteration fallen low enough to stop
    bool changeBelowThreshold = false;
    int numIter = 0;
    // Main loop: update the mesh and check whether it's reached convergence
    while (!changeBelowThreshold) {
        updateMeshGauss(mesh, xSize, ySize, xStep, yStep);

        double maxLoss = 0;
        for (int k = 0; k < numElements; k++) {
            double loss = abs(mesh[k] - oldMesh[k]);
            if (loss > maxLoss)
                maxLoss = loss;
        }

        if (maxLoss < threshold)
            changeBelowThreshold = true; 
        copyMesh(mesh, oldMesh, numElements);
        numIter++;
    }

    if (verbose)
        cout << "Final mesh from rank " << myRank << :\n";

    printMesh(mesh, xSize, ySize);
   
    // Print out extra info about the calculation in verbose mode
    if (verbose) {
        cout << "Expected final mesh:\n";
        printExactResult(xSize, ySize, xStep, yStep);
        // double loss = checkResult(mesh, xSize, ySize, xStep, yStep);
        // cout << "Loss: " << loss << "\n";
        cout << "Number of iterations: " << numIter << "\n";
        // double elapsedTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        // cout << "Total elapsed time = " << elapsedTime << " seconds\n";
    }
}

// Update the values in the mesh using the Gauss-Seidel method
void updateMeshGauss(double* mesh, int xSize, int ySize, double xStep, double yStep) {
    double xStep2 = xStep*xStep;
    double yStep2 = yStep*yStep;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            double source = sourceTerm(((double) j)*xStep, ((double) i)*yStep);
            double term1 = xStep2 * (mesh[(i+1)*xSize + j] + mesh[(i-1)*xSize + j]);
            double term2 = yStep2 * (mesh[i*xSize + (j+1)] + mesh[i*xSize + (j-1)]);
            double term3 = -1 * (xStep2*yStep2*source);
            double result = (term1 + term2 + term3) / (2*(xStep2 + yStep2));
            mesh[i*xSize + j] = result;
        }
    }
}

void updateBounds(double* mesh) {
    myLeftBound = new double[ySize];
    myRightBound = new double[ySize];
    for (int i = 0; i < ySize; i++)
        myLeftBound[i] = mesh[i*
    

// Calculate the source term
double sourceTerm(double x, double y) {
    return x * exp(y);
}

// Print the mesh to standard output
void printMesh(double* mesh, int xSize, int ySize) {
    for (int i = ySize-1; i >= 0; i--) {
        for (int j = 0; j < xSize; j++) 
            cout << setprecision(4) << setw(9) << mesh[i*xSize + j] << " ";
        
        cout << "\n";
    }
    cout << "\n";
}

// Initialize the mesh with 0's and boundary conditions
void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep) {
//            for (int i = 0; i < ySize; i++) {
//            for (int j = 0; j < xSize; j++) {
//                // Bottom bound
//                if (i == 0)
//                    mesh[i*xSize + j] = ((double) j)*xStep;
//                // Top bound
//                else if (i == ySize - 1)
//                    mesh[i*xSize + j] = ((double) j)*xStep*exp(1);
//                // Left bound
//                else if (j == 0)
//                    mesh[i*xSize + j] = 0;
//                // Right bound
//                else if (j == xSize - 1)
//                    mesh[i*xSize + j] = 2*exp(((double) i)*yStep);
//                // Initialize all else to 0
//                else
//                    mesh[i*xSize + j] = 0;
//            }
//        }
//    }
    int numElements = xSize*ySize;
    for (int k = 0; k < numElements; k++) 
        mesh[k] = 0;
}

// Copy one mesh into another
void copyMesh(double* mesh, double* newMesh, int numElements) {
    for (int k = 0; k < numElements; k++) 
        newMesh[k] = mesh[k];
}

// In verify mode, check the result against the exact solution using the max norm
double checkResult(double* mesh, int xSize, int ySize, double xStep, double yStep) {
    double maxLoss = 0;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            double groundTruth = (((double) j)*xStep) * exp(((double) i)*yStep);
            double loss =  abs(groundTruth - mesh[i*xSize + j]);
            if (loss > maxLoss)
                maxLoss = loss;
        }
    }
    return maxLoss;
}

// Print the exact solution
void printExactResult(int xSize, int ySize, double xStep, double yStep) {
    double* exactMesh = new double[xSize*ySize];
    for (int i = 0; i < ySize; i++) {
        for (int j = 0; j < xSize; j++) {
            exactMesh[i*xSize + j] = ((double) j)*xStep * exp(((double) i)*yStep);
        }
    }
    printMesh(exactMesh, xSize, ySize);
    delete exactMesh;
}

void sendBoundary(double* sendBound, int myRank, int xP, int yP, Direction dir) {
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    if (dir == Direction::UP && yPos > 0)
        int destRank = xP*(yPos - 1) + xPos;
    else if (dir == Direction::DOWN && yPos < yP - 1)
        int destRank = xP*(yPos + 1) + xPos;
    else if (dir == Direction::LEFT && xPos > 0)
        int destRank = xP*yPos + xPos - 1;
    else if (dir == Direction::RIGHT && xPos < xP - 1)
        int destRank = xP*yPos + xPos + 1;
    else
        return;
    MPI_Send(sendBound, length(sendBound), MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD);
}

void receiveBoundary(double* recvBound, int myRank, int xP, int yP, int xSize, int ySize, 
                     double xStep, double yStep, Direction dir) {
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    if (dir == Direction::UP) {
        if (yPos == 0) {
            getBoundaryCondition(recvBound, dir, xStep, xSize, xPos)
            return;
        }
        int fromRank = xP*(yPos - 1) + xPos;
    else if (dir == Direction::DOWN)
        if (yPos == yP - 1) {
            getBoundaryCondition(recvBound, dir, xStep, xSize, xPos)
            return;
        }
        int fromRank = xP*(yPos + 1) + xPos;
    else if (dir == Direction::LEFT)
        if (xPos == 0) {
            getBoundaryCondition(recvBound, dir, yStep, ySize, yPos)
            return;
        }
        int fromRank = xP*yPos + xPos - 1;
    else
        if (xPos == xP - 1) {
            getBoundaryCondition(recvBound, dir, yStep, ySize, yPos)
            return;
        }
        int fromRank = xP*yPos + xPos + 1;
    MPI_Recv(&recvBound, length(recvBound), MPI_DOUBLE, fromRank, 0, MPI_COMM_WORLD);
}

void getBoundaryCondition(double* bound, Direction dir, double step, int size, int pos) {
    double val = size*step*pos;
    if (dir == Direction::UP) {
        for (int i = 0; i < size; i++) {
            bound[i] = x*exp(1);
            val += step;
        }
    }
    else if (dir == Direction::DOWN) {
        for (int i = 0; i < size; i++) {
            bound[i] = val;
            val += step;
        }
    }
    else if (dir == Direction::LEFT) {
        for (int i = 0; i < size; i++)
            bound[i] = 0;
    }
    else {
        for (int i = 0; i < size; i++) {
            bound[i] = 2*exp(val);
            val += step;
        }
    }
}
