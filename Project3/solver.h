void updateMeshJacobi(long myRank);
void* work(void* rank);
double sourceTerm(double x, double y);
void printMesh(double* mesh, int xSize, int ySize);
void setBounds(double* mesh);
void copyMesh(double* mesh, double* newMesh, int numElements);
double checkResult();
void printExactResult();
void outputMesh(int myRank, double* mesh, int xSize, int ySize);
