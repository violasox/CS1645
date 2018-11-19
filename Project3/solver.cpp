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
pthread_cond_t condVar;
int xP;
int yP;
int totalXSize;
int totalYSize;
double xStep, yStep;
const double threshold = pow(10, -12);
double* mesh;
double* oldMesh;
double maxLoss;
int numIter = 0;
// Array to monitor the convergence on a solution
double convergence[10000];
int globalCounter = 0;
bool globalStop = false;

// Main function: calculates 2D Poisson with either Jacobi or Gauss-Seidel based on boolean flags and command line arguments
// First two arguments should be size of the mesh in (x y) order
// Optional 3rd argument -v triggers verbose printout mode
int main(int argc, char** argv) {

    if (argc < 5) {
        cout << "Need at least 4 arguments: mesh size in x and y directions and number of processors in x and y directions. Use 5th argument -v for verbose mode.\n";
        return 1;
    }

    totalXSize = stoi(argv[1]);
    totalYSize = stoi(argv[2]);
    xP = stoi(argv[3]);
    yP = stoi(argv[4]);

    if ((totalXSize - 2) % xP != 0) {
        cout << "x mesh size must be evenly divisible by the number of processorsi\n";
        return 1;
    }
    else if ((totalYSize - 2) % yP != 0) {
        cout << "y mesh size must be evenly divisible by the number of processors\n";
        return 1;
    }

    clock_t start = clock();
    // Size of the mesh for each processor
    // Subtract 2 because of the outside boundary
    int xSize = ((totalXSize - 2) / xP);
    int ySize = ((totalYSize - 2) / yP);

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
    
    xStep = ((double) (xMax - xMin)) / ((double) (totalXSize-1));
    yStep = ((double) (yMax - yMin)) / ((double) (totalYSize-1));
    
    if (verbose) 
        cout << "X size: " << xSize << ", Y size: " << ySize << "\n"; 
    
    // Allocate the mesh grid that Poisson should be calculated over
    // Include border around outside for storing neighboring values
    int numElements = totalXSize*totalYSize;
    mesh = new double[numElements];
    setBounds(mesh);
    oldMesh = new double[numElements];
    copyMesh(mesh, oldMesh, numElements);
    printMesh(mesh, totalXSize, totalYSize);
    printMesh(oldMesh, totalXSize, totalYSize);

    maxLoss = 0;
    pthread_mutex_init(&guard, NULL);
    pthread_cond_init(&condVar, NULL);

    int numThreads = xP*yP;
    pthread_t* threadHandles = new pthread_t[numThreads];
    for (long i = 0; i < numThreads; i++)
        pthread_create(&threadHandles[i], NULL, work, (void*) i);

    for (long i = 0; i < numThreads; i++)
        pthread_join(threadHandles[i], NULL);

    double elapsed = (double) (clock() - start) / CLOCKS_PER_SEC;
    if (verbose) {
        cout << "Final mesh:\n";
        printMesh(mesh, totalXSize, totalYSize);
    }
   
    // Print out extra info about the calculation in verbose mode
    if (verbose) {
        cout << "Expected final mesh:\n";
        printExactResult();
        double loss = checkResult();
        cout << "Loss: " << loss << "\n";
        cout << "Number of iterations: " << numIter << "\n";
        cout << "Total elapsed time: " << elapsed << " seconds\n";
    }
    
    // outputMesh(myRank, mesh, xSize, ySize);
    ofstream fout("convergence.txt");
    for (int i = 0; i < numIter / 1000; i++)
        fout << convergence[i] << ",";
}

void* work(void* rank) {
    long myRank = (long) rank;
    int xSize = ((totalXSize - 2) / xP);
    int ySize = ((totalYSize - 2) / yP);
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    void* myReturn;
    // Main loop: update the mesh and check whether it's reached convergence
    while (true) {
        // cout << "Got here " << globalCounter << "\n";
        updateMeshJacobi(myRank);
        //if (numIter % 5 == 0) {
           // cout << "Mesh at iteration " << numIter << " for rank " << myRank << ":\n";
           // printMesh(mesh, xSize, ySize);
        //}
        double myMaxLoss = 0;
        // Calculate the max loss for this thread's section
        for (int i = xSize*xPos; i < xSize*(xPos + 1); i++) {
            for (int j = ySize*yPos; j < ySize*(yPos + 1); j++) {
                double loss = abs(mesh[i*totalXSize + j] - oldMesh[i*totalXSize + j]);
                if (loss > myMaxLoss)
                    myMaxLoss = loss;
            }
        }
        // Barrier function to synchronize threads
        pthread_mutex_lock(&guard);
        if (globalStop) {
            pthread_mutex_unlock(&guard);
            return myReturn;
        }
        if (myMaxLoss > maxLoss)
            maxLoss = myMaxLoss;
        globalCounter++;
        if (globalCounter == (xP*yP)) {
            if (numIter % 1000 == 0)
                convergence[numIter / 1000] = maxLoss;
            globalCounter = 0;
            numIter++;
            copyMesh(mesh, oldMesh, totalXSize*totalYSize);
            if (maxLoss < threshold && numIter > 1) {
                globalStop = true;        
            }
            maxLoss = 0;
            pthread_cond_broadcast(&condVar);
        }
        else {
            while (pthread_cond_wait(&condVar, &guard));
        }
        pthread_mutex_unlock(&guard);  
    }
}


// Update the values in the mesh using the Gauss-Seidel method
void updateMeshJacobi(int myRank) {
    int xSize = ((totalXSize - 2) / xP);
    int ySize = ((totalYSize - 2) / yP);
    double xStep2 = xStep*xStep;
    double yStep2 = yStep*yStep;
    int xPos = myRank % xP;
    int yPos = myRank / xP;
    double xStart = xStep*(xPos*xSize + 1);
    double yStart = yStep*(yPos*ySize + 1);
    for (int i = xSize*xPos + 1; i < xSize*(xPos + 1) + 1; i++) {
        for (int j = ySize*yPos + 1; j < ySize*(yPos + 1) + 1; j++) {
            double source = sourceTerm(((double) j)*xStep + xStart, ((double) i)*yStep + yStart);
            double term1 = xStep2 * (oldMesh[(i+1)*totalXSize + j] + oldMesh[(i-1)*totalXSize + j]);
            double term2 = yStep2 * (oldMesh[i*totalXSize + (j+1)] + oldMesh[i*totalXSize + (j-1)]);
            double term3 = -1 * (xStep2*yStep2*source);
            double result = (term1 + term2 + term3) / (2*(xStep2 + yStep2));
            mesh[i*totalXSize + j] = result;
        }
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
void setBounds(double* mesh) {
    for (int i = 0; i < totalYSize; i++) {
        for (int j = 0; j < totalXSize; j++) {
            // Bottom bound
            if (i == 0)
                mesh[i*totalXSize + j] = ((double) j)*xStep;
            // Top bound
            else if (i == totalYSize - 1)
                mesh[i*totalXSize + j] = ((double) j)*xStep*exp(1);
            // Left bound
            else if (j == 0)
                mesh[i*totalXSize + j] = 0;
            // Right bound
            else if (j == totalXSize - 1)
                mesh[i*totalXSize + j] = 2*exp(((double) i)*yStep);
            // Initialize all else to 0
            else
                mesh[i*totalXSize + j] = 0;
        }
    }
}

// Copy first mesh into the second
void copyMesh(double* mesh, double* newMesh, int numElements) {
    for (int k = 0; k < numElements; k++) 
        newMesh[k] = mesh[k];
}

// In verify mode, check the result against the exact solution using the max norm
double checkResult() {
    double maxErr = 0;
    for (int i = 1; i < totalYSize-1; i++) {
        for (int j = 1; j < totalXSize-1; j++) {
            double groundTruth = (((double) j)*xStep) * exp(((double) i)*yStep);
            double err = abs(groundTruth - mesh[i*totalXSize + j]);
            if (err > maxErr)
                maxErr = err;
        }
    }
    return maxErr;
}

// Print the exact solution
void printExactResult() {
    double* exactMesh = new double[totalXSize*totalYSize];
    for (int i = 0; i < totalYSize; i++) {
        for (int j = 0; j < totalXSize; j++) {
            exactMesh[i*totalXSize + j] = ((double) j)*xStep * exp(((double) i)*yStep);
        }
    }
    printMesh(exactMesh, totalXSize, totalYSize);
    delete exactMesh;
}

