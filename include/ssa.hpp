#ifndef SSA_H
#define SSA_H

#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <vector>

#include "gpgmp.hpp"

void ssa(
        species& x,
        double t_end,
        const std::vector<propensity>& propensities,
        const std::vector<std::vector<int>>& stoich_matrix,
        std::mt19937& gen);

#endif
