double* updateMeshJacobi(double* mesh, int xSize, int ySize);
void updateMeshGauss(double* mesh, double* newMesh, int xSize, int ySize, double xStep, double yStep);
double sourceTerm(int x, int y);
void printMesh(double* mesh, int xSize, int ySize);
void copyMesh(double* mesh, double* newMesh, int numElements);
