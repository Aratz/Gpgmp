#include <QtTest>
#include <QObject>

#include<iostream>

#include <cmath>
#include <cuda.h>

#include "multiparticle.cuh"

class MultiparticleGpuTest: public QObject {
    Q_OBJECT

        private slots :
            void testEmpty();
            void testConstant();
            void testTwoSpecies();
            void testFixed();
            void testUniform();
};

void MultiparticleGpuTest::testEmpty() {
    int Nx = 3;
    int Ny = 3;
    int Ns = 1;
    int shape[4] = {Nx, Ny, Ns, 0};

    grid domain(boost::extents[Nx][Ny][Ns]);

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
    cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    __constant__  int dev_domain_shape[4];
    cudaMemcpyToSymbol(dev_domain_shape, shape, 4 * sizeof(int));

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    for (auto i = 0; i < 100; ++i)
        multiparticle(dev_domain, dev_lost_particles, 0.5, 0, dev_domain_shape);

    cudaMemcpy(&domain[0][0][0], dev_domain, domain.num_elements() * sizeof(int),
            cudaMemcpyDeviceToHost);
    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            QCOMPARE(domain[i][j][0], 0);

    cudaFree(dev_domain);
    cudaFree(dev_lost_particles);
}

void MultiparticleGpuTest::testConstant() {
    int Nx = 3;
    int Ny = 3;
    int Ns = 1;
    int shape[4] = {Nx, Ny, Ns, 0};

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            domain[i][j][0] = 1;

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
    cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    __constant__  int dev_domain_shape[4];
    cudaMemcpyToSymbol(dev_domain_shape, shape, 4 * sizeof(int));

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    for (auto i = 0; i < 100; ++i)
        multiparticle(dev_domain, dev_lost_particles, 0.5, 0, dev_domain_shape);

    cudaMemcpy(&domain[0][0][0], dev_domain, domain.num_elements() * sizeof(int),
            cudaMemcpyDeviceToHost);
    auto sum = 0;
    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            sum += domain[i][j][0];

    QCOMPARE(sum, Nx * Ny);

    cudaFree(dev_domain);
    cudaFree(dev_lost_particles);
}

void MultiparticleGpuTest::testTwoSpecies() {
    int Nx = 3;
    int Ny = 3;
    int Ns = 2;
    int shape[4] = {Nx, Ny, Ns, 0};

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (grid::index s=0; s != Ns; s++)
        for (grid::index i = 0; i != Nx; ++i)
            for (grid::index j = 0; j != Ny; ++j)
                domain[i][j][s] = s + 1;

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
    cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    __constant__  int dev_domain_shape[4];
    cudaMemcpyToSymbol(dev_domain_shape, shape, 4 * sizeof(int));

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    for (auto i = 0; i < 100; ++i)
        for (auto s = 0; s < Ns; ++s)
            multiparticle(dev_domain, dev_lost_particles, 0.5, s, dev_domain_shape);

    cudaMemcpy(&domain[0][0][0], dev_domain, domain.num_elements() * sizeof(int),
            cudaMemcpyDeviceToHost);

    for (grid::index s=0; s != Ns; s++){
        auto sum = 0;
        for (grid::index i = 0; i != Nx; ++i)
            for (grid::index j = 0; j != Ny; ++j)
                sum += domain[i][j][s];
        QCOMPARE(sum, Nx * Ny * ((int) s+1));
    }

    cudaFree(dev_domain);
    cudaFree(dev_lost_particles);
}

void MultiparticleGpuTest::testFixed() {
    int Nx = 3;
    int Ny = 3;
    int Ns = 1;
    int shape[4] = {Nx, Ny, Ns, 0};

    grid domain(boost::extents[Nx][Ny][Ns]);

    domain[0][0][0] = 1;

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
    cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    __constant__  int dev_domain_shape[4];
    cudaMemcpyToSymbol(dev_domain_shape, shape, 4 * sizeof(int));

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    for (auto i = 0; i < 100; ++i)
        multiparticle(dev_domain, dev_lost_particles, 0.5, 0, dev_domain_shape);

    cudaMemcpy(&domain[0][0][0], dev_domain, domain.num_elements() * sizeof(int),
            cudaMemcpyDeviceToHost);

    for (grid::index i = 0; i != Nx; ++i){
        for (grid::index j = 0; j != Ny; ++j){
            if (i || j)
                QCOMPARE(domain[i][j][0], 0);
            else
                QCOMPARE(domain[i][j][0], 1);
        }
    }

    cudaFree(dev_domain);
    cudaFree(dev_lost_particles);
}

void MultiparticleGpuTest::testUniform() {
    int Nx = 3;
    int Ny = 3;
    int Ns = 1;
    int shape[4] = {Nx, Ny, Ns, 0};

    grid domain(boost::extents[Nx][Ny][Ns]);

    int total = 9000;
    domain[0][0][0] = total;

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
    cudaMemcpy(dev_domain, &domain[0][0][0], domain.num_elements() * sizeof(int),
            cudaMemcpyHostToDevice);

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    for (auto i = 0; i < 100; ++i)
        multiparticle(dev_domain, dev_lost_particles, 0.5, 0, shape);

    cudaMemcpy(&domain[0][0][0], dev_domain, domain.num_elements() * sizeof(int),
            cudaMemcpyDeviceToHost);
    for (grid::index i = 0; i != Nx; ++i){
        for (grid::index j = 0; j != Ny; ++j){
            QVERIFY(std::abs(domain[i][j][0] - total/(Nx*Ny)) < total/(Nx*Ny)/5.);
        }
    }


    cudaFree(dev_domain);
    cudaFree(dev_lost_particles);
}

QTEST_MAIN(MultiparticleGpuTest)
#include "MultiparticleGpuTest.moc"
