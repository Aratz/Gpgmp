#include "ssa.hpp"

void ssa(
        species& x,
        double t_end,
        const std::vector<propensity>& propensities,
        const std::vector<species>& stoich_matrix,
        std::mt19937& gen) {
    double t = 0.0;

    while (t < t_end) {
        std::vector<double> a;
        std::uniform_real_distribution<> dis(0.0, 1.0);

        a.reserve(propensities.size());

        std::transform(propensities.begin(), propensities.end(), back_inserter(a),
                [&x] (auto p) { return p(x); });

        std::partial_sum(a.begin(), a.end(), a.begin());

        double a_0 = a.back();

        double tau = log(1./(1.-dis(gen)))/a_0;

        if ((t += tau) < t_end) {
            int i = std::distance(
                    a.begin(),
                    std::lower_bound(a.begin(), a.end(), dis(gen) * a_0));

            std::transform(x.begin(), x.end(), stoich_matrix[i].cbegin(),
                    x.begin(), std::plus<>{});
        }
    }
}
