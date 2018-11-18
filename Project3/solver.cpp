#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <string>
#include <ctime>
#include "solver.h"

using namespace std;
pthread_mutex_t guard;
int xP;
int yP;
int xSize;
int ySize;

// Main function: calculates 2D Poisson with either Jacobi or Gauss-Seidel based on boolean flags and command line arguments
// First two arguments should be size of the mesh in (x y) order
// Optional 3rd argument -v triggers verbose printout mode
int main(int argc, char** argv) {
    double startTime = MPI_Wtime();

    if (argc < 5) {
        cout << "Need at least 4 arguments: mesh size in x and y directions and number of processors in x and y directions. Use 5th argument -v for verbose mode.\n";
        return 1;
    }

    int totalXSize = stoi(argv[1]);
    int totalYSize = stoi(argv[2]);
    xP = stoi(argv[3]);
    yP = stoi(argv[4]);

    if ((totalXSize - 2) % xP != 0) {
        if (myRank == 0)
            cout << "x mesh size must be evenly divisible by the number of processorsi\n";
        return 1;
    }
    else if ((totalYSize - 2) % yP != 0) {
        if (myRank == 0)
            cout << "y mesh size must be evenly divisible by the number of processors\n";
        return 1;
    }
    else if (numP != xP * yP) {
        if (myRank == 0)
            cout << "Total # processors must equal x-direction processors * y-direction processors\n";
        return 1;
    }

    // Size of the mesh for each processor (not including boudnaries)
    // Subtract 2 because of the outside boundary
    int xSize = ((totalXSize - 2) / xP) + 2;
    int ySize = ((totalYSize - 2) / yP) + 2;

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
    if (verbose && myRank == 0)
        cout << "X size: " << xSize << ", Y size: " << ySize << "\n"; 
    
    // Allocate the mesh grid that Poisson should be calculated over
    // Include border around outside for storing neighboring values
    int numElements = xSize*ySize;
    double* mesh = new double[numElements];
    for (int k = 0; k < numElements; k++)
        mesh[k] = 0;
    double* oldMesh = new double[numElements];
    copyMesh(mesh, oldMesh, numElements);

    // Array to monitor the convergence on a solution
    double convergence[10000];

    // Has the change after each iteration fallen low enough to stop
    bool changeBelowThreshold = false;
    bool loop = true;
    bool sentFinished = false;
    bool stopSignal = false;
    int numIter = 0;
    int stopIter = 0;
    MPI_Request broadRequest;
    // Main loop: update the mesh and check whether it's reached convergence
    while (loop) {
        updateMeshGauss(mesh, myRank, xSize, ySize, xP, yP, xStep, yStep);
        //if (numIter % 5 == 0) {
           // cout << "Mesh at iteration " << numIter << " for rank " << myRank << ":\n";
           // printMesh(mesh, xSize, ySize);
        //}
        
        double maxLoss = 0;
        for (int k = 0; k < numElements; k++) {
            double loss = abs(mesh[k] - oldMesh[k]);
            if (loss > maxLoss)
                maxLoss = loss;
        }
        if (myRank == 0 && numIter % 1000 == 0)
            convergence[numIter / 1000] = maxLoss;

        if (maxLoss < threshold) {
            if (stopSignal) {
                if (stopIter == numIter)
                    loop = false;
            }
            else if (myRank != 0) {
                if (!sentFinished) {
                    MPI_Isend(&numIter, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &sendRequest);
                    MPI_Ibcast(&stopIter, 1, MPI_INT, 0, MPI_COMM_WORLD, &broadRequest);
                    sentFinished = true;
                }
                int flag = 0;
                MPI_Test(&broadRequest, &flag, MPI_STATUS_IGNORE);
                if (flag) {
                    stopSignal = true;
                    cout << "Process " << myRank << " received stop signal\n";
                }
            }
            else {
                bool allFinished = true;
                int flag = 0;
                for (int i = 1; i < numP; i++) {
                    MPI_Iprobe(i, 1, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
                    if (flag == 0) {
                        allFinished = false;
                        break;
                    }
                }
                if (allFinished && !sentFinished) {
                    cout << "All processes have finished!\n";
                    stopIter = numIter + 6;
                    cout << "Stop after iteration " << stopIter << "\n";
                    MPI_Ibcast(&stopIter, 1, MPI_INT, 0, MPI_COMM_WORLD, &broadRequest);
                    sentFinished = true;
                    stopSignal = true;
                }
            }
        }
        copyMesh(mesh, oldMesh, numElements);
        numIter++;
        // cout << "Process " << myRank << " finished iteration " << numIter << "\n";
    }
    if (verbose)
        cout << "Final mesh from rank " << myRank << ":\n";

    // printMesh(mesh, xSize, ySize);
   
    // Print out extra info about the calculation in verbose mode
    if (verbose) {
        if (myRank == 0) {
            //cout << "Expected final mesh:\n";
            //printExactResult(totalXSize, totalYSize, xStep, yStep);
        }
        double loss = checkResult(mesh, myRank, xP, xSize, ySize, xStep, yStep);
        cout << "Loss: " << loss << "\n";
        cout << "Number of iterations for rank " << myRank << ": " << numIter << "\n";
        double elapsedTime = MPI_Wtime() - startTime;
        cout << "Total elapsed time for process " << myRank << " is " << elapsedTime << " seconds\n";
    }
    
    outputMesh(myRank, mesh, xSize, ySize);
    if (myRank == 0) {
        ofstream fout("convergence.txt");
        for (int i = 0; i < numIter / 1000; i++)
            fout << convergence[i] << ",";
    }
    MPI_Finalize();
}

// Update the values in the mesh using the Gauss-Seidel method
void updateMeshGauss(double* mesh, int myRank, int xSize, int ySize, int xP, int yP, 
                     double xStep, double yStep) {
    double xStep2 = xStep*xStep;
    double yStep2 = yStep*yStep;
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    double xStart = xStep*xPos*(xSize-2);
    double yStart = yStep*yPos*(ySize-2);
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            double source = sourceTerm(((double) j)*xStep + xStart, ((double) i)*yStep + yStart);
            double term1 = xStep2 * (mesh[(i+1)*xSize + j] + mesh[(i-1)*xSize + j]);
            double term2 = yStep2 * (mesh[i*xSize + (j+1)] + mesh[i*xSize + (j-1)]);
            double term3 = -1 * (xStep2*yStep2*source);
            double result = (term1 + term2 + term3) / (2*(xStep2 + yStep2));
            mesh[i*xSize + j] = result;
        }
    }
    updateBounds(mesh, myRank, xSize, ySize, xP, yP, xStep, yStep);
}

void updateBounds(double* mesh, int myRank, int xSize, int ySize, int xP, int yP,
                  double xStep, double yStep) {
    // Share this processor's updated values (non-blocking)
    double myLeftBound[ySize];
    double myRightBound[ySize];
    for (int i = 0; i < ySize; i++) {
        myLeftBound[i] = mesh[i*xSize + 1];
        myRightBound[i] = mesh[(i+1)*(xSize) - 2];
    }
    sendBoundary(myLeftBound, ySize, myRank, xP, yP, Direction::LEFT);
    sendBoundary(myRightBound, ySize, myRank, xP, yP, Direction::RIGHT);

    double myTopBound[xSize];
    double myBotBound[xSize];
    for (int j = 0; j < xSize; j++) {
        myTopBound[j] = mesh[xSize + j];
        myBotBound[j] = mesh[(ySize-2)*xSize + j];
    }
    sendBoundary(myTopBound, xSize, myRank, xP, yP, Direction::UP);
    sendBoundary(myBotBound, xSize, myRank, xP, yP, Direction::DOWN);
    
    // Receive updated values from the other processors (blocking)
    double otherLeftBound[ySize];
    double otherRightBound[ySize];
    receiveBoundary(otherLeftBound, myRank, xP, yP, xSize, ySize, xStep, yStep, Direction::LEFT);
    receiveBoundary(otherRightBound, myRank, xP, yP, xSize, ySize, xStep, yStep, Direction::RIGHT);
    //cout << "Right boundary: ";
    for (int i = 0; i < ySize; i++) {
        //cout << otherRightBound[i] << ": " << (i+1)*(xSize) - 1 << ", ";
        mesh[i*xSize] = otherLeftBound[i];
        mesh[(i+1)*(xSize) - 1] = otherRightBound[i];
    }
    //cout << "\n";

    double otherTopBound[xSize];
    double otherBotBound[xSize];
    receiveBoundary(otherTopBound, myRank, xP, yP, xSize, ySize, xStep, yStep, Direction::UP);
    receiveBoundary(otherBotBound, myRank, xP, yP, xSize, ySize, xStep, yStep, Direction::DOWN);
    for (int j = 0; j < xSize; j++) {
        mesh[j] = otherTopBound[j];
        mesh[(ySize-1)*xSize + j] = otherBotBound[j];
    }
}
    

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

// Print the mesh to a file
void outputMesh(int myRank, double* mesh, int xSize, int ySize) {
    string filename = to_string(myRank) + "_output.txt";
    ofstream fout(filename);
    for (int i = ySize-2; i >= 1; i--) {
        for (int j = 1; j < xSize-1; j++) 
            fout << setprecision(4) << setw(9) << mesh[i*xSize + j] << " ";    
        fout << "\n";
    }
    fout.close();
}


// Initialize the mesh with 0's, including boundaries
void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep) {
    int numElements = (xSize+2)*(ySize+2);
    for (int k = 0; k < numElements; k++) 
        mesh[k] = 0;
}

// Copy one mesh into another
void copyMesh(double* mesh, double* newMesh, int numElements) {
    for (int k = 0; k < numElements; k++) 
        newMesh[k] = mesh[k];
}

// In verify mode, check the result against the exact solution using the max norm
double checkResult(double* mesh, int myRank, int xP, int xSize, int ySize, double xStep, double yStep) {
    double maxLoss = 0;
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    double xVal;
    double yVal = (double) (yPos*(ySize-2) + 1)*yStep;
    for (int i = 1; i < ySize-1; i++) {
        xVal = (double) (xPos*(xSize-2) + 1)*xStep;
        for (int j = 1; j < xSize-1; j++) {
            double groundTruth = xVal * exp(yVal);
            double loss = abs(groundTruth - mesh[i*xSize + j]);
            if (loss > maxLoss)
                maxLoss = loss;
            xVal += xStep;
        }
        yVal += yStep;
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

void sendBoundary(double* sendBound, int boundLen, int myRank, int xP, int yP, Direction dir) {
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    int destRank;
    if (dir == Direction::UP && yPos > 0)
        destRank = xP*(yPos - 1) + xPos;
    else if (dir == Direction::DOWN && yPos < yP - 1)
        destRank = xP*(yPos + 1) + xPos;
    else if (dir == Direction::LEFT && xPos > 0)
        destRank = xP*yPos + xPos - 1;
    else if (dir == Direction::RIGHT && xPos < xP - 1)
        destRank = xP*yPos + xPos + 1;
    else
        return;
    //cout << "Sending boundary\n";
    MPI_Isend(sendBound, boundLen, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &sendRequest);
}

void receiveBoundary(double* recvBound, int myRank, int xP, int yP, int xSize, int ySize, 
                     double xStep, double yStep, Direction dir) {
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    int boundLen, fromRank;
    if (dir == Direction::UP) {
        if (yPos == 0) {
            getBoundaryCondition(recvBound, dir, xStep, xSize, xPos);
            return;
        }
        fromRank = xP*(yPos - 1) + xPos;
        boundLen = xSize;
    }
    else if (dir == Direction::DOWN) {
        if (yPos == yP - 1) {
            getBoundaryCondition(recvBound, dir, xStep, xSize, xPos);
            return;
        }
        fromRank = xP*(yPos + 1) + xPos;
        boundLen = xSize;
    }
    else if (dir == Direction::LEFT) {
        if (xPos == 0) {
            getBoundaryCondition(recvBound, dir, yStep, ySize, yPos);
            return;
        }
        fromRank = xP*yPos + xPos - 1;
        boundLen = ySize;
    }
    else {
        if (xPos == xP - 1) {
            getBoundaryCondition(recvBound, dir, yStep, ySize, yPos);
            return;
        }
        fromRank = xP*yPos + xPos + 1;
        boundLen = ySize;
    }
    MPI_Recv(recvBound, boundLen, MPI_DOUBLE, fromRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// If getting left/right boundary, use yStep, ySize, yPos
// Opposite if getting up/down boundary
void getBoundaryCondition(double* bound, Direction dir, double step, int size, int pos) {
    // cout << "Receiving boundary condition for ";
    // Find the initial x/y value for this processor's section
    // Subtract 2 because repeated unit doesn't include boundaries
    double val = (size-2)*step*pos;
    if (dir == Direction::DOWN) {
        //cout << "up\n";
        for (int i = 0; i < size; i++) {
            bound[i] = val*exp(1);
            val += step;
        }
    }
    else if (dir == Direction::UP) {
        //cout << "down\n";
        for (int i = 0; i < size; i++) {
            bound[i] = val;
            val += step;
        }
    }
    else if (dir == Direction::LEFT) {
        //cout << "left\n";
        for (int i = 0; i < size; i++)
            bound[i] = 0;
    }
    else {
        //cout << "right\n";
        for (int i = 0; i < size; i++) {
            bound[i] = 2*exp(val);
            val += step;
        }
    }
}
