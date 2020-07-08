#ifndef GPGMP_H
#define GPGMP_H

#include <algorithm>

#include "grid.hpp"
#include "ssa.hpp"
#include "multiparticle.hpp"

void simulate(
        grid& domain,
        double t_end,
        const std::vector<propensity>& propensities,
        const std::vector<std::vector<int>>& stoich_matrix,
        double p,
        const std::vector<double>& d_s,
        double lambda,
        std::mt19937& gen
        );

#endif
