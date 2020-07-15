#include <QtTest>
#include <QObject>

#include <cmath>

#include "multiparticle.cuh"

class MultiparticleGpuTest: public QObject {
    Q_OBJECT

        private slots :
            void testEmpty();
            /*void testConstant();
            void testTwoSpecies();
            void testFixed();
            void testUniform();
            //*/
};

void MultiparticleGpuTest::testEmpty() {
    int Nx = 3;
    int Ny = 3;
    int Ns = 1;
    int shape[4] = {Nx, Ny, Ns, 0};

    int* dev_domain;
    cudaMalloc((void **) &dev_domain, Nx * Ny * Ns * sizeof(int));
    for (auto i = 0; i < Nx * Ny * Ns; ++i)
        dev_domain[i] = 0;

    __constant__  int dev_domain_shape[4];
    cudaMemcpyToSymbol(dev_domain_shape, shape, 4 * sizeof(int));

    int* dev_lost_particles;
    cudaMalloc((void **) &dev_lost_particles, 4 * Nx * Ny * sizeof(int));

    for (auto i = 0; i < 100; ++i)
        multiparticle(dev_domain, dev_lost_particles, 0.5, 0, dev_domain_shape);

    for (auto i = 0; i < Nx * Ny * Ns; ++i)
        QCOMPARE(dev_domain[i],  0);
}

/*
void MultiparticleGpuTest::testConstant() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            domain[i][j][0] = 1;

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.5, 0, gen);

    auto sum = 0;
    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            sum += domain[i][j][0];

    QCOMPARE(sum, Nx * Ny);
}

void MultiparticleGpuTest::testTwoSpecies() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 2;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (grid::index s=0; s != Ns; s++)
        for (grid::index i = 0; i != Nx; ++i)
            for (grid::index j = 0; j != Ny; ++j)
                domain[i][j][s] = s + 1;

    for (auto i = 0; i < 100; ++i)
        for (grid::index s=0; s != Ns; s++)
            multiparticle(domain, 0.5, s, gen);

    for (grid::index s=0; s != Ns; s++){
        auto sum = 0;
        for (grid::index i = 0; i != Nx; ++i)
            for (grid::index j = 0; j != Ny; ++j)
                sum += domain[i][j][s];
        QCOMPARE(sum, Nx * Ny * ((int) s+1));
    }
}

void MultiparticleGpuTest::testFixed() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    domain[0][0][0] = 1;

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.0, 0, gen);

    for (grid::index i = 0; i != Nx; ++i){
        for (grid::index j = 0; j != Ny; ++j){
            if (i || j)
                QCOMPARE(domain[i][j][0], 0);
            else
                QCOMPARE(domain[i][j][0], 1);
        }
    }
}

void MultiparticleGpuTest::testUniform() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    int total = 9000;
    domain[0][0][0] = total;

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.5, 0, gen);

    for (grid::index i = 0; i != Nx; ++i){
        for (grid::index j = 0; j != Ny; ++j){
            QVERIFY(std::abs(domain[i][j][0] - total/(Nx*Ny)) < total/(Nx*Ny)/10.);
        }
    }
}
//*/
QTEST_MAIN(MultiparticleGpuTest)
#include "MultiparticleGpuTest.moc"
