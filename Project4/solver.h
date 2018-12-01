void updateMeshJacobi(double* mesh, double* newMesh, int xSize, int ySize, double xStep, double yStep, bool verifyMode);
void updateMeshGauss(double* mesh, int xSize, int ySize, double xStep, double yStep, bool verifyMode);
double sourceTerm(double x, double y, bool verifyMode);
void printMesh(double* mesh, int xSize, int ySize);
void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep, bool verifyMode);
void copyMesh(double* mesh, double* newMesh, int numElements);
double checkResult(double* mesh, int xSize, int ySize, double xStep, double yStep);
void printExactResult(int xSize, int ySize, double xStep, double yStep);
