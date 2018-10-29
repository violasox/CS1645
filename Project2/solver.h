enum class Direction { UP, DOWN, LEFT, RIGHT };

void updateMeshGauss(double* mesh, int myRank, int xSize, int ySize, int xP, int yP, 
                     double xStep, double yStep);
double sourceTerm(double x, double y);
void printMesh(double* mesh, int xSize, int ySize);
void setBounds(double* mesh, int xSize, int ySize, double xStep, double yStep, 
               bool verifyMode);
void copyMesh(double* mesh, double* newMesh, int numElements);
double checkResult(double* mesh, int xSize, int ySize, double xStep, double yStep);
void printExactResult(int xSize, int ySize, double xStep, double yStep);
void sendBoundary(double* sendBound, int myRank, int xP, int yP, Direction dir);
void receiveBoundary(double* recvBound, int myRank, int xP, int yP, int xSize,
                     int ySize, double xStep, double yStep, Direction dir);
void getBoundaryCondition(double* bound, Direction dir, double step, int size,
                          int pos);
void updateBounds(double* mesh, int myRank, int xSize, int ySize, int xP, int yP,
                  double xStep, double yStep);

