#include <QtTest>
#include <QObject>

#include <cmath>
#include <random>
#include <vector>

#include "gpgmp.hpp"

class GpgmpTest: public QObject {
    Q_OBJECT

        private slots :
            void testGpgmp();
            void testTime();
};

void GpgmpTest::testGpgmp() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 2;

    grid domain(boost::extents[Nx][Ny][Ns]);
    int n = 900;
    domain[1][1][0] = n;

    double t_end = 1000.;

    double mu = 1000.;
    const std::vector<propensity> propensities {
        [mu] (const species& x) { return mu*x[0]; },
    };

    std::vector<std::vector<int>> stoich_matrix {
        { -1, 1 }
    };

    const std::vector<double> d_s = { 1., 1. };
    double lambda = 1.;
    double p = 0.5;

    simulate(
            domain,
            t_end,
            propensities,
            stoich_matrix,
            p,
            d_s,
            lambda,
            gen);

    int a = 0, b = 0;
    for (auto x = 0; x < Nx; ++x) {
        for (auto y = 0; y < Ny; ++y) {
            QVERIFY(std::abs(domain[x][y][1] - 1.*n/(Nx*Ny)) < 0.25*n/(Nx*Ny));
            a += domain[x][y][0];
            b += domain[x][y][1];
        }
    }

    QCOMPARE(a, 0);
    QCOMPARE(b, n);
}

void GpgmpTest::testTime() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int Nx = 3;
    int Ny = 3;
    int Ns = 1;

    grid domain(boost::extents[Nx][Ny][Ns]);

    double t_end = 500.;

    const std::vector<double> d_s = { 0.001 };
    double lambda = 0.5;
    double p = 0.5;


    double mu = 1./(Nx*Ny);
    const std::vector<propensity> propensities {
        [mu] (const species& x) { return mu; },
    };

    std::vector<std::vector<int>> stoich_matrix {
        { 1 }
    };

    simulate(
            domain,
            t_end,
            propensities,
            stoich_matrix,
            p,
            d_s,
            lambda,
            gen);

    int a = 0;
    for (auto x = 0; x < Nx; ++x)
        for (auto y = 0; y < Ny; ++y)
            a += domain[x][y][0];

    QVERIFY(std::abs(a - t_end) < 0.2*t_end);
}

QTEST_MAIN(GpgmpTest)
#include "GpgmpTest.moc"
