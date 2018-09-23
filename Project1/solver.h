double* updateMeshJacobi(double* mesh, int xSize, int ySize);
void updateMeshGauss(double* mesh, double* newMesh, int xSize, int ySize, double xStep, double yStep);
double sourceTerm(double x, double y);
void printMesh(double* mesh, int xSize, int ySize);
void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep);
void copyMesh(double* mesh, double* newMesh, int numElements);
