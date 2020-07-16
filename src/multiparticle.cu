#include "multiparticle.cuh"

__device__ int * subdomain(int* domain, int x, int y, int s, int * domain_shape){
    int Ny = domain_shape[1];
    int Ns = domain_shape[2];

    return domain + Ny*Ns*x + Ns*y + s;
}

__device__ int * sublost_particles(int* lost_particles, int x, int y, int i, int * domain_shape){
    int Ny = domain_shape[1];

    return lost_particles + 4*Ny*x + 4*y + i;
}


__global__ void multiparticlegpu(
        int* domain,
        int* lost_particles,
        double p,
        int s,
        int* domain_shape
        ) {
    int threadID;
    threadID = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int seed = threadID;
    curandState crs;
    curand_init(seed, 0, 0, &crs);

    int Nx = domain_shape[0];
    int Ny = domain_shape[1];

    const int x = threadIdx.x;
    const int y = threadIdx.y;

    for(int i = 0; i < 4; ++i){
        *sublost_particles(lost_particles, x, y, i, domain_shape) = 0;
    }

    int dn = 0;
    for (int n = 0; n < *subdomain(domain, x, y, s, domain_shape); ++n){
        int u = curand_uniform(&crs) < p;
        *sublost_particles(lost_particles,
                x, y, (int) (curand_uniform(&crs)*2*DIM),
                domain_shape) += u;
        dn += u;
    }

    *subdomain(domain, x, y, s, domain_shape) -= dn;

    __syncthreads();

    for (int d = 0; d < 2*DIM; d++) {
        int dx = (d % DIM == 0) * (d < DIM?-1:1);
        int dy = (d % DIM == 1) * (d < DIM?-1:1);

        for (int x = 0; x < Nx; ++x){
            for (int y = 0; y < Ny; ++y){
                int outofbounds = (0 > x + dx || Nx <= x + dx
                                    || 0 > y + dy || Ny <= y + dy);

                if (!outofbounds) {
                    *subdomain(domain, x, y, s, domain_shape) += *sublost_particles(
                            lost_particles, x + dx, y + dy, d, domain_shape);
                }
                else {
                    *subdomain(domain, x, y, s, domain_shape) += *sublost_particles(
                            lost_particles, 
                            (x < -dx || x + dx >= Nx?x:x+dx),
                            (y < -dy || y + dy >= Ny?y:y+dy),
                            (d + DIM)%(2*DIM),
                            domain_shape);
                }
            }
        }
    }
}

__host__ void multiparticle(
        int* domain,
        int* lost_particles,
        double p,
        int s,
        int domain_shape[4]
        ) {
    int Nx = domain_shape[0];
    int Ny = domain_shape[1];


    int* dev_domain_shape;
    cudaMalloc((void **) &dev_domain_shape, 4 * sizeof(int));
    cudaMemcpy(dev_domain_shape, domain_shape, 4 * sizeof(int),
            cudaMemcpyHostToDevice);

    dim3 N(Nx, Ny, 1);

    multiparticlegpu<<<1, N>>>(domain, lost_particles, 0.5, 0, dev_domain_shape);

    cudaFree(dev_domain_shape);
}
