int main(int xSize, int ySize, double xMin, double xMax, double yMin, double yMax) {
    // Allocate the mesh grid that Poisson should be calculated over
    double** mesh = new double*[xSize];
    for (int i = 0; i < xSize; i++) {
        mesh[i] = new int[ySize];
        // Set the initial conditions
        for (int j = 0; j < ySize; j++) {
            mesh[i][j] = 0;
            if (i == 0)
                mesh[i][j] = xMin;
            else if (i == xSize - 1)
                mesh[i][j] = xMax;
            else if (j == 0)
                mesh[i][j] = yMin;
            else if (j == ySize - 1)
                mesh[i][j] = yMax;
    }

    // Has the change after each iteration fallen low enough to stop
    bool changeBelowThreshold = false;
    // Test value to check change
    double oldTestVal = 0;
    const double threshold = 0.001
    while (!changeBelowThreshold) {
        updateMeshJacobi(mesh, double xMin, double xMax, double yMin, double yMax)
        newTestVal = mesh[xSize/2][ySize/2];
        if (abs(newTestVal - oldTestVal) < threshold)
            changeBelowThreshold = true; 
    }

double** updateMeshJacobi(mesh, double xMin, double xMax, double yMin, double yMax) {
    return mesh;
}

double** updateMeshGauss(mesh, double xMin, double xMax, double yMin, double yMax) {
    return mesh;
}
