#include "ssa.hpp"

__device__ int* get_species(int* domain, int x, int y, int s, int* domain_shape){
    int Nx = domain_shape[0];
    int Ny = domain_shape[1];
    int Ns = domain_shape[2];

    return domain[Ny*Ns*x + Ns*y + s];
}

__device__ int* get_rates(int* texture_array, int x , int y, int r, int* domain_shape) {
    int Nx = domain_shape[0];
    int Ny = domain_shape[1];
    int Nr = domain_shape[3];

    return texture_array + Ny*Nr*x + Nr*y + r;

}

__device__ double propensity(int* domain, int x, int y, int r,
        int* texture_array, int* domain_shape) {
    int* rates = get_rates(texture_array, x, y, r, domain_shape);
    int w = rates[0];
    int k = rates[1];
    int Ni = rates[2];
    int Nj = rates[3];

    return w * k
        * (Ni >= 0?*get_species(domain, x, y, Ni, domain_shape):1)
        * (Nj >= 0?*get_species(domain, x, y, Nj, domain_shape):1);
}

__global__ void ssa(
        int* domain,
        int* domain_shape,
        double t_end,
        int* texture_array,
        int* stoich_matrix
        ) {

    int Nx = domain_shape[0];
    int Ny = domain_shape[1];
    int Ns = domain_shape[2];
    int Nr = domain_shape[3];

    if(Nr == 0){
        return;
    }

    const int x = blockIdx.x;
    const int y = blockIdx.y;

    int threadID;
    threadID = blockIdx.y * blockDim.x + threadIdx.x;
    unsigned int seed = threadID + t_end*1000;
    curandState crs;
    curand_init(seed, 0, 0, &crs);

    double t = 0.0;

    while (t < *t_end) {
        double* a;
        cudaMalloc((void **) &a, Nr * sizeof(double));
        a[0] = propensity(domain, x, y, 0, domain_shape);
        for (int i = 1; i < Nr; ++i){
            a[i] = a[i-1] + propensity(domain, x, y, i, domain_shape);
        }

        double a_0 = a[Nr - 1];

        double tau = log(1./(1.-curand_uniform(&crs)))/a_0;

        if ((t += tau) < t_end) {
            int r = 0;
            double u = curand_uniform(&crs) * a_0;

            while(r + 1 < Nr && a[r + 1] < u)
                r += 1;

            for (int s = 0; s < Ns; s++) {
                *get_species(domain, x, y, s, domain_shape) += stoich_matrix[Ns*r + s];
            }
        }
    }

    cudaFree(a);
}
