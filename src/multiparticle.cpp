#include "multiparticle.hpp"

void multiparticle(
        grid& domain,
        double p,
        int s,
        std::mt19937& gen) {
    std::uniform_real_distribution<> dis(0.0, 1.0);

    auto shape = domain.shape();
    auto Nx = shape[0];
    auto Ny = shape[1];
    auto Ns = shape[2];

    boost::multi_array<int, DIM + 2> lost_particles(boost::extents[Nx][Ny][Ns][2*DIM]);

    for (auto x = 0; x < Nx; ++x){
        for (auto y = 0; y < Ny; ++y){
            int dn = 0;
            for (auto n = 0; n < domain[x][y][s]; ++n){
                int u = dis(gen) < p;
                lost_particles[x][y][s][(size_t) (dis(gen)*2*DIM)] += u;
                dn += u;
            }
            domain[x][y][s] -= dn;
        }
    }

    for (auto d = 0; d < 2*DIM; d++) {
        int dx = (d % DIM == 0) * (d < DIM?-1:1);
        int dy = (d % DIM == 1) * (d < DIM?-1:1);

        for (auto x = 0; x < Nx; ++x){
            for (auto y = 0; y < Ny; ++y){
                auto outofbounds = (0 > x + dx || Nx <= x + dx
                                    || 0 > y + dy || Ny <= y + dy);

                if (!outofbounds) {
                    domain[x][y][s] += lost_particles[x + dx][y + dy][s][d];
                }
                else {
                    /**/
                    domain[x][y][s] += lost_particles
                        [(x < -dx || x + dx >= Nx?x:x+dx)]
                        [(y < -dy || y + dy >= Ny?y:y+dy)]
                        [s]
                        [(d + DIM)%(2*DIM)];
                }
            }
        }
    }
}
