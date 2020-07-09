#include "ssa.hpp"

void ssa(
        grid& domain,
        double t_end,
        const std::vector<propensity>& propensities,
        const std::vector<std::vector<int>>& stoich_matrix,
        std::mt19937& gen) {

    auto shape = domain.shape();
    auto Nx = shape[0];
    auto Ny = shape[1];
    auto Ns = shape[2];

    for (auto x = 0; x < Nx; ++x){
        for (auto y = 0; y < Ny; ++y){
            double t = 0.0;
            auto subvolume = domain[boost::indices[x][y][range(0, Ns)]];

            while (t < t_end) {
                std::vector<double> a;
                std::uniform_real_distribution<> dis(0.0, 1.0);

                a.reserve(propensities.size());

                std::transform(propensities.begin(), propensities.end(), back_inserter(a),
                        [&subvolume] (auto p) { return p(subvolume); });

                std::partial_sum(a.begin(), a.end(), a.begin());

                double a_0 = a.back();

                double tau = log(1./(1.-dis(gen)))/a_0;

                if ((t += tau) < t_end) {
                    int i = std::distance(
                            a.begin(),
                            std::lower_bound(a.begin(), a.end(), dis(gen) * a_0));

                    std::transform(subvolume.begin(), subvolume.end(), stoich_matrix[i].cbegin(),
                            subvolume.begin(), std::plus<>{});
                }
            }
        }
    }
}
