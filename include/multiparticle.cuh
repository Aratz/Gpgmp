#ifndef CUMP_H
#define CUMP_H

#include <cuda.h>
#include <curand_kernel.h>

#include "grid.hpp"

__global__ void multiparticle(
        int* domain,
        int* lost_particles,
        double p,
        int s,
        int* domain_shape);

#endif
