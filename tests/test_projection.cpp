#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/operators.h"
#include "ns2d/projection.h"
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
