#include <iostream>
#include <vector>
#include <chrono>

#include <cuda.h>

#include "multiparticle.cuh"

int main() {
    std::vector<int> gridsizes  {10, 10, 100, 1000, 10000};

    for (int gridsize: gridsizes) {
        int Nx = gridsize;
        int Ny = gridsize;
        int Ns = 1;
        int shape[4] = {Nx, Ny, Ns, 0};

        grid domain(boost::extents[Nx][Ny][Ns]);

        for (int x = 0; x < gridsize; ++x)
            for (int y = 0; y < gridsize; ++y)
                domain[0][0][0] = 1000;

        auto begin = std::chrono::steady_clock::now();

        int* dev_domain;
        cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
        cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
                cudaMemcpyHostToDevice);

        int* dev_lost_particles;
        cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

        for (auto i = 0; i < 10000; ++i)
            multiparticle(dev_domain, dev_lost_particles, 0.5, 0, shape);

        cudaMemcpy(&domain[0][0][0], dev_domain, domain.num_elements() * sizeof(int),
                cudaMemcpyDeviceToHost);

        cudaFree(dev_domain);
        cudaFree(dev_lost_particles);

        auto end = std::chrono::steady_clock::now();

        std::cout << "Benchmark " << gridsize << "x" << gridsize << ": "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
            << " microseconds\n";
    }
}
