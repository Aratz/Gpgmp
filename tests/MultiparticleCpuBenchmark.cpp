#include <iostream>
#include <QtTest>
#include <QObject>

#include <cmath>
#include <random>
#include <vector>

#include "multiparticle.hpp"


class MultiparticleCpuBenchmark: public QObject {
    Q_OBJECT

        private slots :
            void benchmarkMultiparticle64();
            void benchmarkMultiparticle128();
            void benchmarkMultiparticle256();
};

void MultiparticleCpuBenchmark::benchmarkMultiparticle64() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int gridsize = 64;

    int Nx = gridsize;
    int Ny = gridsize;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (int x = 0; x < gridsize; ++x)
        for (int y = 0; y < gridsize; ++y)
            domain[x][y][0] = 10;

    QBENCHMARK {
        for (auto i = 0; i < 100; ++i)
            multiparticle(domain, 0.5, 0, gen);
    }
}

void MultiparticleCpuBenchmark::benchmarkMultiparticle128() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int gridsize = 128;

    int Nx = gridsize;
    int Ny = gridsize;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (int x = 0; x < gridsize; ++x)
        for (int y = 0; y < gridsize; ++y)
            domain[x][y][0] = 10;

    QBENCHMARK {
        for (auto i = 0; i < 100; ++i)
            multiparticle(domain, 0.5, 0, gen);
    }
}

void MultiparticleCpuBenchmark::benchmarkMultiparticle256() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int gridsize = 256;

    int Nx = gridsize;
    int Ny = gridsize;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (int x = 0; x < gridsize; ++x)
        for (int y = 0; y < gridsize; ++y)
            domain[x][y][0] = 10;

    QBENCHMARK {
        for (auto i = 0; i < 100; ++i)
            multiparticle(domain, 0.5, 0, gen);
    }
}

QTEST_MAIN(MultiparticleCpuBenchmark)
#include "MultiparticleCpuBenchmark.moc"
