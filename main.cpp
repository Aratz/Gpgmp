#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

typedef std::vector<int> species;
typedef std::function<double(const species&)> propensity;

double ssa(
        species& x,
        const std::vector<propensity>& propensities,
        const std::vector<species>& stoich_matrix,
        std::mt19937& gen) {
    std::vector<double> a;
    std::uniform_real_distribution<> dis(0.0, 1.0);

    a.reserve(propensities.size());

    std::transform(propensities.cbegin(), propensities.cend(), back_inserter(a),
            [&x] (auto p) { return p(x); });

    std::partial_sum(a.begin(), a.end(), a.begin());

    double a_0 = a.back();
    double tau = log(1./(1.-dis(gen)))/a_0;

    int i = std::distance(
            a.cbegin(),
            std::lower_bound(a.cbegin(), a.cend(), dis(gen) * a_0));

    std::transform(x.cbegin(), x.cend(), stoich_matrix[i].cbegin(),
            x.begin(), std::plus<>{});

    return tau;
}

int main() {
    return 0;
}
