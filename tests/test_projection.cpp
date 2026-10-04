#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/operators.h"
#include "ns2d/projection.h"
#include "ns2d/poisson.h"
#include "ns2d/boundary.h"
#include <algorithm>
#include <random>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <numbers>

namespace {

constexpr double pi = std::numbers::pi;

// Fill every stored point (ghosts included) with func at the u or v location
template <typename Func>
void fillU(Field& f, const Grid& g, Func func) {
    for (int j = 0; j <= g.getNj() - 1; ++j) {
        for (int i = 0; i <= g.getNi() - 1; ++i) {
            f(i, j) = func(g.xu(i), g.yp(j));
        }
    }
}

template <typename Func>
void fillV(Field& f, const Grid& g, Func func) {
    for (int j = 0; j <= g.getNj() - 1; ++j) {
        for (int i = 0; i <= g.getNi() - 1; ++i) {
            f(i, j) = func(g.xp(i), g.yv(j));
        }
    }
}

bool insideRange(const IndexRange& r, int i, int j) {
    return i >= r.iBegin && i <= r.iEnd && j >= r.jBegin && j <= r.jEnd;
}

} // namespace

TEST_CASE("Predictor keeps fluid at rest", "[projection]") {
    const Grid g(6, 4, 2.0, 1.0);
    const Field u(g, 0.0);
    const Field v(g, 0.0);
    Field uStar(g, -999.0);
    Field vStar(g, -999.0);
    PredictorWork work(g);

    predictor(u, v, g, 0.01, 100.0, AdvectionScheme::Upwind2, work, uStar, vStar);

    // Exact comparison: nothing is computed from non-zero values
    CHECK(uStar.values() == u.values());
    CHECK(vStar.values() == v.values());
}

TEST_CASE("Predictor leaves uniform flow unchanged", "[projection]") {
    const Grid g(6, 4, 2.0, 1.0);
    const Field u(g, 2.5);
    const Field v(g, 0.0);
    Field uStar(g);
    Field vStar(g);
    PredictorWork work(g);

    for (AdvectionScheme scheme : {AdvectionScheme::Central, AdvectionScheme::Upwind2}) {
        INFO("scheme = " << static_cast<int>(scheme));
        predictor(u, v, g, 0.01, 100.0, scheme, work, uStar, vStar);

        // Advection and Laplacian of a constant field are 0 up to rounding
        for (std::size_t k = 0; k < u.size(); ++k) {
            CHECK_THAT(uStar.values()[k], Catch::Matchers::WithinAbs(2.5, 1e-12));
            CHECK_THAT(vStar.values()[k], Catch::Matchers::WithinAbs(0.0, 1e-12));
        }
    }
}

TEST_CASE("Predictor combines operators with the right signs and factors", "[projection]") {
    const Grid g(16, 16, 2.0, 1.0);
    const double dt = 0.01;
    const double Re = 50.0;

    Field u(g);
    Field v(g);
    fillU(u, g, [](double x, double y) { return std::sin(pi * x) * std::cos(pi * y); });
    fillV(v, g, [](double x, double y) { return std::cos(pi * x) * std::sin(2.0 * pi * y); });

    Field uStar(g);
    Field vStar(g);
    PredictorWork work(g);
    predictor(u, v, g, dt, Re, AdvectionScheme::Upwind2, work, uStar, vStar);

    // Expected: compute the operators independently of the predictor
    Field advU(g), advV(g), lapU(g), lapV(g);
    advectionU(u, v, g, AdvectionScheme::Upwind2, advU);
    advectionV(u, v, g, AdvectionScheme::Upwind2, advV);
    laplacian(u, g, Staggering::UFace, lapU);
    laplacian(v, g, Staggering::VFace, lapV);

    const IndexRange ru = interiorRange(g, Staggering::UFace);
    const IndexRange rv = interiorRange(g, Staggering::VFace);

    SECTION("interior points follow u* = u + dt * (-adv + lap / Re)") {
        for (int j = ru.jBegin; j <= ru.jEnd; ++j) {
            for (int i = ru.iBegin; i <= ru.iEnd; ++i) {
                INFO("u: i = " << i << ", j = " << j);
                const double expected = dt * (-advU(i, j) + lapU(i, j) / Re);
                CHECK_THAT(uStar(i, j) - u(i, j), Catch::Matchers::WithinAbs(expected, 1e-12));
            }
        }
        for (int j = rv.jBegin; j <= rv.jEnd; ++j) {
            for (int i = rv.iBegin; i <= rv.iEnd; ++i) {
                INFO("v: i = " << i << ", j = " << j);
                const double expected = dt * (-advV(i, j) + lapV(i, j) / Re);
                CHECK_THAT(vStar(i, j) - v(i, j), Catch::Matchers::WithinAbs(expected, 1e-12));
            }
        }
    }

    SECTION("walls and ghosts keep the old values") {
        for (int j = 0; j <= g.getNj() - 1; ++j) {
            for (int i = 0; i <= g.getNi() - 1; ++i) {
                INFO("i = " << i << ", j = " << j);
                if (!insideRange(ru, i, j)) {
                    CHECK(uStar(i, j) == u(i, j));
                }
                if (!insideRange(rv, i, j)) {
                    CHECK(vStar(i, j) == v(i, j));
                }
            }
        }
    }
}

namespace {

double maxAbsInterior(const Field& f, const Grid& g) {
    double m = 0.0;
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            m = std::max(m, std::abs(f(i, j)));
        }
    }
    return m;
}

} // namespace

TEST_CASE("Divergence is exact for linear velocity fields", "[projection]") {
    const Grid g(6, 4, 2.0, 1.0);
    Field u(g);
    Field v(g);
    fillU(u, g, [](double x, double) { return 3.0 * x; });
    fillV(v, g, [](double, double y) { return -2.0 * y + 1.0; });

    Field div(g, -999.0);
    divergence(u, v, g, div);

    // du/dx + dv/dy = 3 - 2 = 1 at every cell
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            INFO("i = " << i << ", j = " << j);
            CHECK_THAT(div(i, j), Catch::Matchers::WithinAbs(1.0, 1e-12));
        }
    }
    // Ghosts of the output are not written
    CHECK(div(0, 1) == -999.0);
    CHECK(div(g.getNx() + 1, 1) == -999.0);
}

TEST_CASE("Projection removes the divergence of a random velocity field", "[projection]") {
    const Grid g(16, 12, 2.0, 1.0);   // non-square, dx != dy
    const int nx = g.getNx();
    const int ny = g.getNy();
    const double dt = 0.01;

    // Random u*, v*, then closed walls: normal wall velocities become 0
    Field uStar(g);
    Field vStar(g);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    std::generate(uStar.values().begin(), uStar.values().end(), [&]() { return dis(gen); });
    std::generate(vStar.values().begin(), vStar.values().end(), [&]() { return dis(gen); });
    const BoundaryConditions walls{};
    walls.applyNoSlipWall(uStar, vStar, g);

    Field div(g);
    divergence(uStar, vStar, g, div);
    const double divBefore = maxAbsInterior(div, g);

    Field rhs(g, 0.0);
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            rhs(i, j) = div(i, j) / dt;
        }
    }

    Field p(g, 0.0);
    JacobiSolver jacobi(g, 500000);
    const SolveResult res = jacobi.solve(rhs, p, 1e-10);
    INFO("Jacobi iterations = " << res.iterations << ", residual = " << res.residual);
    REQUIRE(res.residual < 1e-10);

    // Start u, v from garbage: the corrector must produce a complete result on its own
    Field u(g, 123.0);
    Field v(g, 123.0);
    corrector(uStar, vStar, p, g, dt, u, v);

    SECTION("divergence is removed down to solver tolerance") {
        divergence(u, v, g, div);
        const double divAfter = maxAbsInterior(div, g);
        INFO("max|div u*| = " << divBefore << ", max|div u| = " << divAfter);
        CHECK(divBefore > 1.0);                   // the test field really was divergent
        CHECK(divAfter / divBefore < 1e-8);
    }

    SECTION("wall normal velocities stay zero") {
        for (int j = 1; j <= ny; ++j) {
            INFO("j = " << j);
            CHECK(u(0, j) == 0.0);
            CHECK(u(nx, j) == 0.0);
        }
        for (int i = 1; i <= nx; ++i) {
            INFO("i = " << i);
            CHECK(v(i, 0) == 0.0);
            CHECK(v(i, ny) == 0.0);
        }
    }

    SECTION("ghosts and walls are copied from u*, v*") {
        const IndexRange ru = interiorRange(g, Staggering::UFace);
        const IndexRange rv = interiorRange(g, Staggering::VFace);
        for (int j = 0; j <= g.getNj() - 1; ++j) {
            for (int i = 0; i <= g.getNi() - 1; ++i) {
                INFO("i = " << i << ", j = " << j);
                if (!insideRange(ru, i, j)) {
                    CHECK(u(i, j) == uStar(i, j));
                }
                if (!insideRange(rv, i, j)) {
                    CHECK(v(i, j) == vStar(i, j));
                }
            }
        }
    }
}

namespace {

// Set every interior u and v value (ghosts are left as they are)
void setInterior(Field& u, Field& v, const Grid& g, double uVal, double vVal) {
    const IndexRange ru = interiorRange(g, Staggering::UFace);
    const IndexRange rv = interiorRange(g, Staggering::VFace);
    for (int j = ru.jBegin; j <= ru.jEnd; ++j) {
        for (int i = ru.iBegin; i <= ru.iEnd; ++i) {
            u(i, j) = uVal;
        }
    }
    for (int j = rv.jBegin; j <= rv.jEnd; ++j) {
        for (int i = rv.iBegin; i <= rv.iEnd; ++i) {
            v(i, j) = vVal;
        }
    }
}

} // namespace

TEST_CASE("computeDt picks the stricter of the advective and diffusive limits", "[projection]") {
    const Grid g(8, 8, 2.0, 1.0);   // dx = 0.25, dy = 0.125
    const double dx = g.getDx();
    const double dy = g.getDy();
    const double cfl = 0.5;
    Field u(g, 0.0);
    Field v(g, 0.0);

    SECTION("fluid at rest: diffusive limit") {
        const double Re = 100.0;
        const double dtDiff = 0.5 * Re / (1.0 / (dx * dx) + 1.0 / (dy * dy));
        CHECK_THAT(computeDt(u, v, g, Re, cfl), Catch::Matchers::WithinRel(dtDiff, 1e-12));
    }

    SECTION("fast x-flow: advective limit") {
        const double Re = 1000.0;
        setInterior(u, v, g, 100.0, 0.0);
        CHECK_THAT(computeDt(u, v, g, Re, cfl), Catch::Matchers::WithinRel(cfl * dx / 100.0, 1e-12));
    }

    SECTION("both components count") {
        const double Re = 1000.0;
        setInterior(u, v, g, 100.0, 50.0);
        const double expected = cfl / (100.0 / dx + 50.0 / dy);
        CHECK_THAT(computeDt(u, v, g, Re, cfl), Catch::Matchers::WithinRel(expected, 1e-12));
    }

    SECTION("negative velocities count by magnitude") {
        const double Re = 1000.0;
        setInterior(u, v, g, -100.0, -50.0);
        const double expected = cfl / (100.0 / dx + 50.0 / dy);
        CHECK_THAT(computeDt(u, v, g, Re, cfl), Catch::Matchers::WithinRel(expected, 1e-12));
    }

    SECTION("doubling the speed halves dt in the advective regime") {
        const double Re = 1000.0;
        setInterior(u, v, g, 100.0, 50.0);
        const double dt1 = computeDt(u, v, g, Re, cfl);
        setInterior(u, v, g, 200.0, 100.0);
        const double dt2 = computeDt(u, v, g, Re, cfl);
        CHECK_THAT(dt2, Catch::Matchers::WithinRel(0.5 * dt1, 1e-12));
    }

    SECTION("ghost values are ignored") {
        const double Re = 1000.0;
        setInterior(u, v, g, 1.0, 0.0);
        const double dtClean = computeDt(u, v, g, Re, cfl);
        u(1, 0) = 1000.0;                       // bottom ghost
        u(1, g.getNy() + 1) = 1000.0;           // lid ghost
        v(0, 1) = 1000.0;                       // left ghost
        CHECK(computeDt(u, v, g, Re, cfl) == dtClean);
    }
}
