#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include <cassert>
#include<cmath>

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
    double mu = 100.0;
    double gamma = 1.0;

    double t_end = 1000.;

    std::random_device rd;
    std::mt19937 gen(rd());

    species x {0};
    std::vector<species> stoich_matrix {
        { 1 },
        { -1 }
    };
    const std::vector<propensity> propensities {
        [mu] (const species& x) { return mu; },
        [gamma] (const species& x) { return gamma * x[0]; },
    };

    double t = 0.0;
    while (t < t_end) {
        std::cout << "t:" << t << "; ";
        std::cout << "P:" << x[0] << "\n";
        t += ssa(x, propensities, stoich_matrix, gen);
    }
    std::cout << "P:" << x[0] << "\n";
    std::cout << "t:" << t << "\n";

    assert(std::abs((x[0] - mu/gamma)/(mu/gamma)) < 0.2);
    return 0;
}
