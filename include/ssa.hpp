#ifndef SSA_H
#define SSA_H

#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <vector>

typedef std::vector<int> species;
typedef std::function<double(const species&)> propensity;

void ssa(
        species& x,
        double t_end,
        const std::vector<propensity>& propensities,
        const std::vector<species>& stoich_matrix,
        std::mt19937& gen);

#endif
