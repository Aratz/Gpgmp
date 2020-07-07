#include <QtTest>
#include <QObject>

#include <cmath>
#include <random>
#include <vector>

#include "multiparticle.hpp"

class MultiparticleTest: public QObject {
    Q_OBJECT

        private slots :
            void testEmpty();
            void testConstant();
            void testTwoSpecies();
            void testFixed();
            void testUniform();
};

void MultiparticleTest::testEmpty() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.5, gen);

    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            QCOMPARE(domain[i][j][0], 0);
}

void MultiparticleTest::testConstant() {
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
        multiparticle(domain, 0.5, gen);

    auto sum = 0;
    for (grid::index i = 0; i != Nx; ++i)
        for (grid::index j = 0; j != Ny; ++j)
            sum += domain[i][j][0];

    QCOMPARE(sum, Nx * Ny);
}

void MultiparticleTest::testTwoSpecies() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 2;

    grid domain(boost::extents[Nx][Ny][Ns]);

    for(grid::index s=0; s!= Ns; s++)
        for (grid::index i = 0; i != Nx; ++i)
            for (grid::index j = 0; j != Ny; ++j)
                domain[i][j][s] = s + 1;

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.5, gen);

    for(grid::index s=0; s!= Ns; s++){
        auto sum = 0;
        for (grid::index i = 0; i != Nx; ++i)
            for (grid::index j = 0; j != Ny; ++j)
                sum += domain[i][j][s];
        QCOMPARE(sum, Nx * Ny * ((int) s+1));
    }
}

void MultiparticleTest::testFixed() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    domain[0][0][0] = 1;

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.0, gen);

    for (grid::index i = 0; i != Nx; ++i){
        for (grid::index j = 0; j != Ny; ++j){
            if (i || j)
                QCOMPARE(domain[i][j][0], 0);
            else
                QCOMPARE(domain[i][j][0], 1);
        }
    }
}

void MultiparticleTest::testUniform() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    int total = 9000;
    domain[0][0][0] = total;

    for (auto i = 0; i < 100; ++i)
        multiparticle(domain, 0.5, gen);

    for (grid::index i = 0; i != Nx; ++i){
        for (grid::index j = 0; j != Ny; ++j){
            QVERIFY(std::abs(domain[i][j][0] - total/(Nx*Ny)) < total/(Nx*Ny)/10.);
        }
    }
}

QTEST_MAIN(MultiparticleTest)
#include "MultiparticleTest.moc"
