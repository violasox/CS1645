#include <cmath>
#include <fstream>
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
    // Flag to use Jacobi (true) or Gauss-Seidel (false)
    bool useJacobi = false;
    
    if (argc < 3) {
        cout << "Need at least arguments: mesh size in x and y directions. Use 3rd argument -v for verbose mode.\n";
        return 1;
    }

    int xSize = stoi(argv[1]);
    int ySize = stoi(argv[2]);

    // Flag to trigger printing out useful info (set in command line with -v)
    bool verbose = false;
    if (argc > 3 && string(argv[3]).compare("-v") == 0)
        verbose = true;

    double xMin, xMax, yMin, yMax;
    if (verifyMode) {
        xMin = 0;
        xMax = 2;
        yMin = 0;
        yMax = 1;
    }
    else {
        xMin = 0;
        xMax = 1;
        yMin = 0;
        yMax = 1;
    }
    
    double xStep = ((double) (xMax - xMin)) / ((double) (xSize-1));
    double yStep = ((double) (yMax - yMin)) / ((double) (ySize-1));
    
    const double threshold = pow(10, -12);
    
    // Allocate the mesh grid that Poisson should be calculated over
    int numElements = xSize*ySize;
    double* mesh = new double[numElements];
    setBounds(mesh, xSize, ySize, xStep, yStep, verifyMode);
    double* mesh2 = new double[numElements];
    copyMesh(mesh, mesh2, numElements);
    double* oldMesh = new double[numElements];
    copyMesh(mesh, oldMesh, numElements);

    if (verbose) {
        //cout << "Starting mesh:\n";
        //printMesh(mesh, xSize, ySize);
    }

    double convergence[10000];

    // Has the change after each iteration fallen low enough to stop
    bool changeBelowThreshold = false;
    int numIter = 0;
    // Main loop: update the mesh and check whether it's reached convergence
    while (!changeBelowThreshold) {
        if (useJacobi)
            updateMeshJacobi(mesh, mesh2, xSize, ySize, xStep, yStep, verifyMode);
        else
            updateMeshGauss(mesh, xSize, ySize, xStep, yStep, verifyMode);

        double maxLoss = 0;
        for (int k = 0; k < numElements; k++) {
            double loss = abs(mesh[k] - oldMesh[k]);
            if (loss > maxLoss)
                maxLoss = loss;
        }
        
        if (numIter % 1000 == 0)
            convergence[numIter / 1000] = maxLoss;

        if (maxLoss < threshold)
            changeBelowThreshold = true; 
        copyMesh(mesh, oldMesh, numElements);
        numIter++;
    }

    if (verbose)
        cout << "Final mesh:\n";

    // printMesh(mesh, xSize, ySize);
   
    // Print out extra info about the calculation in verbose mode
    if (verbose) {
        if (verifyMode) {
            //cout << "Expected final mesh:\n";
            //printExactResult(xSize, ySize, xStep, yStep);
            double loss = checkResult(mesh, xSize, ySize, xStep, yStep);
            cout << "Loss: " << loss << "\n";
        }
        else {
            double testVal = mesh[(ySize/2)*xSize + (xSize/2)];
            cout << "Temperature in middle: " << testVal << "\n";
        }
        cout << "Number of iterations: " << numIter << "\n";
        double elapsedTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        cout << "Total elapsed time = " << elapsedTime << " seconds\n";
    }
    ofstream fout("convergence.txt");
    for (int i = 0; i < numIter / 1000; i++)
        fout << convergence[i] << ",";
}

// Update the values in the mesh using the Gauss-Seidel method
void updateMeshGauss(double* mesh, int xSize, int ySize, double xStep, double yStep, bool verifyMode) {
    double xStep2 = xStep*xStep;
    double yStep2 = yStep*yStep;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            double source = sourceTerm(((double) j)*xStep, ((double) i)*yStep, verifyMode);
            double term1 = xStep2 * (mesh[(i+1)*xSize + j] + mesh[(i-1)*xSize + j]);
            double term2 = yStep2 * (mesh[i*xSize + (j+1)] + mesh[i*xSize + (j-1)]);
            double term3 = -1 * (xStep2*yStep2*source);
            double result = (term1 + term2 + term3) / (2*(xStep2 + yStep2));
            mesh[i*xSize + j] = result;
        }
    }
}

// Update the values in the mesh using the Jacobi method
void updateMeshJacobi(double* mesh, double* newMesh, int xSize, int ySize, double xStep, double yStep, bool verifyMode) {
    double xStep2 = xStep*xStep;
    double yStep2 = yStep*yStep;
    double source, term1, term2, term3, result;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            source = sourceTerm(((double) j)*xStep, ((double) i)*yStep, verifyMode);
            term1 = xStep2 * (mesh[(i+1)*xSize + j] + mesh[(i-1)*xSize + j]);
            term2 = yStep2 * (mesh[i*xSize + (j+1)] + mesh[i*xSize + (j-1)]);
            term3 = -1 * (xStep2*yStep2*source);
            result = (term1 + term2 + term3) / (2*(xStep2 + yStep2));
            newMesh[i*xSize + j] = result;
        }
    }
    copyMesh(newMesh, mesh, xSize*ySize);
}

// Calculate the source term
double sourceTerm(double x, double y, bool verifyMode) {
    if (verifyMode)
        return x * exp(y);
    else
        return 0.2;
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
void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep, bool verifyMode) {
    if (verifyMode) {
        for (int i = 0; i < ySize; i++) {
            for (int j = 0; j < xSize; j++) {
                // Bottom bound
                if (i == 0)
                    mesh[i*xSize + j] = ((double) j)*xStep;
                // Top bound
                else if (i == ySize - 1)
                    mesh[i*xSize + j] = ((double) j)*xStep*exp(1);
                // Left bound
                else if (j == 0)
                    mesh[i*xSize + j] = 0;
                // Right bound
                else if (j == xSize - 1)
                    mesh[i*xSize + j] = 2*exp(((double) i)*yStep);
                // Initialize all else to 0
                else
                    mesh[i*xSize + j] = 0;
            }
        }
    }
    else {
        for (int k = 0; k < xSize*ySize; k++) 
            mesh[k] = 0;
    }
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
