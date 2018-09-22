#include <cmath>
#include <iostream>
#include <string>
#include "solver.h"

using namespace std;

int main(int argc, char** argv) {
    int xSize = 5;
    int ySize = 5;
    int xMin = 0;
    int xMax = 2;
    int yMin = 0;
    int yMax = 1;
    double topBound = 0;
    double botBound = 2;
    double leftBound = 1;
    double rightBound = 3; 
    
    double xStep = ((double) (xMax - xMin)) / ((double) xSize);
    xStep *= xStep;
    double yStep = ((double) (yMax - yMin)) / ((double) ySize);
    yStep *= yStep;
    
    const double threshold = pow(10, -12);
    // Allocate the mesh grid that Poisson should be calculated over
    int numElements = xSize*ySize;
    double* mesh = new double[numElements];
    for (int i = 0; i < numElements; i++) {
        mesh[i] = 0;
        // Set the initial conditions
        if (i < xSize)
            mesh[i] = topBound;
        else if (i >= (ySize-1)*xSize)
            mesh[i] = botBound;
        else if (i % ySize == 0)
            mesh[i] = leftBound;
        else if ((i+1) % ySize == 0)
            mesh[i] = rightBound;
    }

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
    while (!changeBelowThreshold) {
        updateMeshGauss(mesh, mesh2, xSize, ySize, xStep, yStep);
        newTestVal = mesh[testIndex];
        if (abs(newTestVal - oldTestVal) < threshold)
            changeBelowThreshold = true; 
        oldTestVal = newTestVal;
        // printMesh(mesh, xSize, ySize);
    }
    cout << "Final mesh:\n";
    printMesh(mesh, xSize, ySize);
}

double* updateMeshJacobi(double* mesh, int xSize, int ySize) {
    return mesh;
}

void updateMeshGauss(double* mesh, double* newMesh, int xSize, int ySize, double xStep, double yStep) {
    double source, term1, term2, term3, result;
    for (int i = 1; i < ySize-1; i++) {
        for (int j = 1; j < xSize-1; j++) {
            source = sourceTerm(j, i);
            term1 = xStep * (mesh[(i+1)*xSize + j] + mesh[(i-1)*xSize + j]);
            term2 = yStep * (mesh[i*xSize + (j+1)] + mesh[i*xSize + (j-1)]);
            term3 = -1 * (xStep*yStep*source);
            result = (term1 + term2 + term3) / (2*(xStep + yStep));
            newMesh[i*xSize + j] = result;
        }
    }
    copyMesh(newMesh, mesh, xSize*ySize);
}

double sourceTerm(int x, int y) {
    return ((double) x) * exp((double) y);
}

void printMesh(double* mesh, int xSize, int ySize) {
    for (int i = 0; i < ySize; i++) {
        for (int j = 0; j < xSize; j++) 
            cout << mesh[i*xSize + j] << " ";
        
        cout << "\n";
    }
    cout << "\n";
}

// TODO: only copy boundary values
void copyMesh(double* mesh, double* newMesh, int numElements) {
    for (int i = 0; i < numElements; i++) 
        newMesh[i] = mesh[i];
}
