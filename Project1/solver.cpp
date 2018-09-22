#include <cmath>
#include <iostream>
#include <string>
#include "solver.h"

using namespace std;

int main(int argc, char** argv) {
    int xSize = 4;
    int ySize = 4;
    int xMin = 0;
    int xMax = 2;
    int yMin = 0;
    int yMax = 1;
    double topBound = 0;
    double botBound = 2;
    double leftBound = 1;
    double rightBound = 3; 
    
    const double threshold = 0.001;
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

    cout << "Starting matrix\n";
    printMesh(mesh, xSize, ySize);

    // Has the change after each iteration fallen low enough to stop
    bool changeBelowThreshold = true;
    // Test value to check change
    double oldTestVal = 0;
    double newTestVal = 0;
    // Index of the array to test
    int testIndex = (ySize/2)*xSize + (xSize/2);
    while (!changeBelowThreshold) {
        updateMeshJacobi(mesh, xSize, ySize);
        newTestVal = mesh[testIndex];
        if (abs(newTestVal - oldTestVal) < threshold)
            changeBelowThreshold = true; 
    }
}

double* updateMeshJacobi(double* mesh, int xSize, int ySize) {
    return mesh;
}

double* updateMeshGauss(double* mesh, int xSize, int ySize) {

    return mesh;
}

void printMesh(double* mesh, int xSize, int ySize) {
    for (int i = 0; i < ySize; i++) {
        for (int j = 0; j < xSize; j++) {
            cout << mesh[i*xSize + j] << " ";
        }
        cout << "\n";
    }
}
