#include <iostream>
#include <random>
#include <vector>

#include <cassert>
#include<cmath>

#include "ssa.hpp"

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
