__device__ double calcStep(double x, double deltaX);
__global__ void deviceIntegrate(double* result);
__device__ double function(double x);
