#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include "solver.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        cout << "Need exactly two arguments: mesh size in x and y directions\n";
        return 1;
    }

    int xSize = stoi(argv[1]);
    int ySize = stoi(argv[2]);
    int xMin = 0;
    int xMax = 2;
    int yMin = 0;
    int yMax = 1;
    
    double xStep = ((double) (xMax - xMin)) / ((double) (xSize-1));
    double yStep = ((double) (yMax - yMin)) / ((double) (ySize-1));
    
    const double threshold = pow(10, -12);
    
    // Allocate the mesh grid that Poisson should be calculated over
    int numElements = xSize*ySize;
    double* mesh = new double[numElements];
    setBounds(mesh, xSize, ySize, xStep, yStep);
    double* mesh2 = new double[numElements];
    copyMesh(mesh, mesh2, numElements);

    cout << "Starting mesh:\n";
    printMesh(mesh, xSize, ySize);

    // Has the change after each iteration fallen low enough to stop
    bool changeBelowThreshold = false;
    // Test value to check change
    double oldTestVal = 0;
    double newTestVal = 0;
    // Index of the array to test
    int testIndex = (ySize/2)*xSize + (xSize/2);
    int numIter = 0;
    while (!changeBelowThreshold) {
        updateMeshGauss(mesh, mesh2, xSize, ySize, xStep, yStep);
        newTestVal = mesh[testIndex];
        if (abs(newTestVal - oldTestVal) < threshold)
            changeBelowThreshold = true; 
        oldTestVal = newTestVal;
        numIter++;
    }
    cout << "Final mesh:\n";
    printMesh(mesh, xSize, ySize);
    cout << "Expected final mesh:\n";
    printExactResult(xSize, ySize, xStep, yStep);
    double loss = checkResult(mesh, xSize, ySize, xStep, yStep);
    cout << "Loss: " << loss << "\n";
    cout << "Number of iterations: " << numIter << "\n";
}

double* updateMeshJacobi(double* mesh, int xSize, int ySize) {
    return mesh;
}

void updateMeshGauss(double* mesh, double* newMesh, int xSize, int ySize, double xStep, double yStep) {
    double xStep2 = xStep*xStep;
    double yStep2 = yStep*yStep;
    double source, term1, term2, term3, result;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            source = sourceTerm(((double) j)*xStep, ((double) i)*yStep);
            term1 = xStep2 * (mesh[(i+1)*xSize + j] + mesh[(i-1)*xSize + j]);
            term2 = yStep2 * (mesh[i*xSize + (j+1)] + mesh[i*xSize + (j-1)]);
            term3 = -1 * (xStep2*yStep2*source);
            result = (term1 + term2 + term3) / (2*(xStep2 + yStep2));
            newMesh[i*xSize + j] = result;
        }
    }
    copyMesh(newMesh, mesh, xSize*ySize);
}

double sourceTerm(double x, double y) {
    return x * exp(y);
}

void printMesh(double* mesh, int xSize, int ySize) {
    for (int i = ySize-1; i >= 0; i--) {
        for (int j = 0; j < xSize; j++) 
            cout << setprecision(4) << setw(6) << mesh[i*xSize + j] << " ";
        
        cout << "\n";
    }
    cout << "\n";
}

void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep) {
    for (int i = 0; i < ySize; i++) {
        for (int j = 0; j < xSize; j++) {
            // Bottom bound
            if (i == 0)
                mesh[i*xSize + j] = ((double) j)*xStep;
            // Top bound
            else if (i == ySize - 1)
                // mesh[i*xSize + j] = exp(((double) j)*xStep);
                mesh[i*xSize + j] = ((double) j)*xStep*exp(1);
            // Left bound
            else if (j == 0)
                mesh[i*xSize + j] = 0;
            // Right bound
            else if (j == ySize - 1)
                mesh[i*xSize + j] = 2*exp(((double) i)*yStep);
        }
    }
}

void copyMesh(double* mesh, double* newMesh, int numElements) {
    for (int i = 0; i < numElements; i++) 
        newMesh[i] = mesh[i];
}

double checkResult(double* mesh, int xSize, int ySize, double xStep, double yStep) {
    double loss = 0;
    double groundTruth;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            groundTruth = (((double) j)*xStep) * exp(((double) i)*yStep);
            loss += abs(groundTruth - mesh[i*xSize + j]);
        }
    }
    int numPoints = (xSize-2)*(ySize-2);
    loss /= (double) numPoints;
    return loss;
}

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
