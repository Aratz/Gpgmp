# include "gpgmp.hpp"

void simulate(
        grid& domain,
        double t_end,
        const std::vector<propensity>& propensities,
        const std::vector<std::vector<int>>& stoich_matrix,
        double p,
        const std::vector<double>& d_s,
        double lambda,
        std::mt19937& gen
        ){
    //TODO cover case where d_s = 0
    auto shape = domain.shape();
    auto Nx = shape[0];
    auto Ny = shape[1];
    auto Ns = shape[2];

    std::vector<double> dt_s;
    dt_s.reserve(d_s.size());

    std::transform(d_s.begin(), d_s.end(), back_inserter(dt_s),
            [p, lambda] (double d_s) { return ((1 - p)/(2.*DIM))*(lambda * lambda / d_s); });

    std::vector<int> n_s(d_s.size(), 1);

    std::vector<double> t_s(dt_s);

    double t = 0.0;

    while (t < t_end) {
        auto it_s = std::min_element(t_s.begin(), t_s.end());

        for (auto x = 0; x < Nx; ++x){
            for (auto y = 0; y < Ny; ++y){
                auto subvolume = domain[boost::indices[x][y][range(0, Ns)]];
                ssa(
                        subvolume,
                        *it_s - t,
                        propensities,
                        stoich_matrix,
                        gen);
            }
        }
        t = *it_s;

        multiparticle(
                domain,
                p,
                it_s - t_s.begin(),
                gen);

        *it_s = ++n_s[it_s - t_s.begin()] * dt_s[it_s - t_s.begin()];
    }
}
