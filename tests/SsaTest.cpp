#include <QtTest>
#include <QObject>

#include <cmath>
#include <random>
#include <vector>

#include "ssa.hpp"

class SsaTest: public QObject {
    Q_OBJECT

        private slots :
            void testSsa();
};

//TODO test empty
//TODO test no reaction

void SsaTest::testSsa() {
    double mu = 100.0;
    double gamma = 1.0;

    double t_end = 1000.;

    std::random_device rd;
    std::mt19937 gen(rd());

    grid domain(boost::extents[1][1][1]);
    species x = domain[boost::indices[0][0][range(0, 1)]];

    std::vector<std::vector<int>> stoich_matrix {
        { 1 },
        { -1 }
    };
    const std::vector<propensity> propensities {
        [mu] (const species& x) { return mu; },
        [gamma] (const species& x) { return gamma * x[0]; },
    };

    ssa(x, t_end, propensities, stoich_matrix, gen);

    QVERIFY(std::abs((x[0] - mu/gamma)/(mu/gamma)) < 0.2);
}

QTEST_MAIN(SsaTest)
#include "SsaTest.moc"
