# include "gpgmp.hpp"

//TODO rename file back to .cpp
void simulate(
        grid& domain,
        double t_end,
        const boost::multi_array<int, DIM + 2>& texture_array,
        const std::vector<std::vector<int>>& stoich_matrix,
        double p,
        const std::vector<double>& d_s,
        double lambda,
        std::mt19937& gen
        ){
    //TODO cover case where d_s = 0
    auto shape = domain.shape();
    int Nx = shape[0];
    int Ny = shape[1];
    int Ns = shape[2];

    int Nr = stoich_matrix.size();

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, domain.num_elements() * sizeof(int));
    cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    __constant__  int dev_domain_shape[4];
    cudaMemcpyToSymbol(dev_domain_shape, shape, 3 * sizeof(int));
    cudaMemcpyToSymbol(dev_domain_shape+3, &Nr, sizeof(int));

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    int* dev_texture_array;
    cudaMalloc((void **) &dev_texture_array, 4 * Nx * Ny * Nr * sizeof(double));
    cudaMemcpy(dev_texture_array, &texture_array[0][0][0][0],
            texture_array.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    int* dev_stoich_matrix;
    cudaMalloc((void **) &dev_stoich_matrix, Nr * Ns * sizeof(int));
    for (auto r = stoich_matrix.begin(); r != stoich_matrix.end(); ++r){
        cudaMemcpy(dev_stoich_matrix + (r - stoich_matrix.begin())*Ns,
                r, Ns * sizeof(int));
    }

    std::vector<double> dt_s;
    dt_s.reserve(d_s.size());

    std::transform(d_s.begin(), d_s.end(), back_inserter(dt_s),
            [p, lambda] (double d_s) { return ((1 - p)/(2.*DIM))*(lambda * lambda / d_s); });

    std::vector<int> n_s(d_s.size(), 1);

    std::vector<double> t_s(dt_s);

    double t = 0.0;

    while (t < t_end) {
        auto it_s = std::min_element(t_s.begin(), t_s.end());

        ssa<<<N, 1>>>(
                domain,
                *it_s - t,
                propensities,
                stoich_matrix,
                gen);

        t = *it_s;

        dim3 N(Nx, Ny);
        multiparticle<<<N,1>>>(
                dev_domain,
                dev_lost_particles,
                p,
                it_s - t_s.begin(),
                dev_domain_shape,
                );

        *it_s = ++n_s[it_s - t_s.begin()] * dt_s[it_s - t_s.begin()];
    }

    //TODO copy back the results!

    cudaFree(dev_domain);
    cudaFree(dev_stoich_matrix);
    cudaFree(dev_lost_particles);
    cudaFree(dev_texture_array);
}
