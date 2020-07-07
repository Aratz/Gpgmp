#ifndef GPGMP_H
#define GPGMP_H

#include "boost/multi_array.hpp"

#define DIM 2

typedef boost::multi_array<int, DIM + 1> grid;
typedef boost::multi_array_types::index_range range;
typedef grid::array_view<1>::type species;
typedef std::function<double(const species&)> propensity;

#endif
