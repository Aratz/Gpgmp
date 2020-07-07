#ifndef MP_H
#define MP_H

#include <algorithm>
#include <random>
#include <vector>

#include "grid.hpp"

void multiparticle(
        grid& domain,
        double p,
        std::mt19937& gen
        );

#endif
