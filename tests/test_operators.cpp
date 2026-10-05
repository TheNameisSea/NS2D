#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/operators.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>

namespace {

struct Point {
    double x{};
    double y{};
};

// Physical location of index (i, j) for a field with staggering s
Point position(const Grid& g, Staggering s, int i, int j) {
    switch (s) {
    case Staggering::Cell:
        return Point{g.xp(i), g.yp(j)};
    case Staggering::UFace:
        return Point{g.xu(i), g.yp(j)};
    case Staggering::VFace:
        return Point{g.xp(i), g.yv(j)};
    }
    return Point{};
}

// Fill every stored point (ghosts included) with func evaluated at its position
template <typename Func>
void fillExact(Field& f, const Grid& g, Staggering s, Func func) {
    for (int j = 0; j <= g.getNj() - 1; ++j) {
        for (int i = 0; i <= g.getNi() - 1; ++i) {
            const auto [x, y] = position(g, s, i, j);
            f(i, j) = func(x, y);
        }
    }
}

// Manufactured solution, not symmetric in x and y
constexpr double pi = std::numbers::pi;

double exactF(double x, double y) {
    return std::sin(pi * x) * std::cos(2.0 * pi * y);
}

double exactLaplacian(double x, double y) {
    return -5.0 * pi * pi * exactF(x, y);
}

// Max-norm error of the discrete Laplacian on an N x N grid with dx != dy
double maxLaplacianError(int N, Staggering s) {
    const Grid g(N, N, 2.0, 1.0);
    Field f(g);
    Field out(g);
    fillExact(f, g, s, exactF);

    laplacian(f, g, s, out);

    double err = 0.0;
    const auto [i0, i1, j0, j1] = interiorRange(g, s);
    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            const auto [x, y] = position(g, s, i, j);
            err = std::max(err, std::abs(out(i, j) - exactLaplacian(x, y)));
        }
    }
    return err;
}

constexpr Staggering allStaggerings[] = {Staggering::Cell, Staggering::UFace, Staggering::VFace};

} // namespace

TEST_CASE("interiorRange on a 4x3 grid", "[operators]") {
    const Grid g(4, 3, 2.0, 1.0);

    SECTION("Cell") {
        const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::Cell);
        CHECK(i0 == 1);
        CHECK(i1 == 4);
        CHECK(j0 == 1);
        CHECK(j1 == 3);
    }
    SECTION("UFace") {
        const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::UFace);
        CHECK(i0 == 1);
        CHECK(i1 == 3);
        CHECK(j0 == 1);
        CHECK(j1 == 3);
    }
    SECTION("VFace") {
        const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::VFace);
        CHECK(i0 == 1);
        CHECK(i1 == 4);
        CHECK(j0 == 1);
        CHECK(j1 == 2);
    }
}

TEST_CASE("Laplacian of x^2 + y^2 is exactly 4", "[operators]") {
    const Grid g(4, 3, 2.0, 1.0);
    constexpr double marker = -999.0;

    for (Staggering s : allStaggerings) {
        INFO("staggering = " << static_cast<int>(s));

        // Fresh fields per staggering, so no values leak between iterations
        Field f(g);
        Field out(g, marker);
        fillExact(f, g, s, [](double x, double y) { return x * x + y * y; });

        laplacian(f, g, s, out);

        const auto [i0, i1, j0, j1] = interiorRange(g, s);
        for (int j = j0; j <= j1; ++j) {
            for (int i = i0; i <= i1; ++i) {
                CHECK_THAT(out(i, j), Catch::Matchers::WithinAbs(4.0, 1e-10));
            }
        }

        // Points just outside the range must be untouched
        CHECK(out(i1 + 1, j0) == marker);
        CHECK(out(i0, j1 + 1) == marker);
        CHECK(out(i0 - 1, j0) == marker);
        CHECK(out(i0, j0 - 1) == marker);
    }
}

TEST_CASE("Laplacian converges with second order", "[operators]") {
    for (Staggering s : allStaggerings) {
        INFO("staggering = " << static_cast<int>(s));

        const double e16 = maxLaplacianError(16, s);
        const double e32 = maxLaplacianError(32, s);
        const double e64 = maxLaplacianError(64, s);
        INFO("errors: " << e16 << " " << e32 << " " << e64);

        CHECK_THAT(std::log2(e16 / e32), Catch::Matchers::WithinAbs(2.0, 0.1));
        CHECK_THAT(std::log2(e32 / e64), Catch::Matchers::WithinAbs(2.0, 0.1));
    }
}

TEST_CASE("Interpolation between u and v points is exact for linear fields", "[operators]") {
    const Grid g(4, 3, 2.0, 1.0);

    // Averaging 4 surrounding points is exact for any linear function
    auto linear = [](double x, double y) { return 3.0 * x + 2.0 * y + 1.0; };

    SECTION("v interpolated to u points") {
        Field v(g);
        fillExact(v, g, Staggering::VFace, linear);

        const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::UFace);
        for (int j = j0; j <= j1; ++j) {
            for (int i = i0; i <= i1; ++i) {
                INFO("i = " << i << ", j = " << j);
                const auto [x, y] = position(g, Staggering::UFace, i, j);
                CHECK_THAT(vAtU(v, i, j), Catch::Matchers::WithinAbs(linear(x, y), 1e-12));
            }
        }
    }

    SECTION("u interpolated to v points") {
        Field u(g);
        fillExact(u, g, Staggering::UFace, linear);

        const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::VFace);
        for (int j = j0; j <= j1; ++j) {
            for (int i = i0; i <= i1; ++i) {
                INFO("i = " << i << ", j = " << j);
                const auto [x, y] = position(g, Staggering::VFace, i, j);
                CHECK_THAT(uAtV(u, i, j), Catch::Matchers::WithinAbs(linear(x, y), 1e-12));
            }
        }
    }
}

TEST_CASE("derivative is exact for quadratics", "[operators]") {
    const Grid g(4, 3, 2.0, 1.0);
    const Staggering s = Staggering::UFace;

    // Central and second-order upwind are both exact for quadratics,
    // so any wrong coefficient, sign, branch or spacing shows up here.
    // Wrong fallback logic shows up as an index assert in Debug builds.
    Field fx(g);
    Field fy(g);
    fillExact(fx, g, s, [](double x, double) { return x * x; });
    fillExact(fy, g, s, [](double, double y) { return y * y; });

    const auto [i0, i1, j0, j1] = interiorRange(g, s);

    for (AdvectionScheme scheme : {AdvectionScheme::Central, AdvectionScheme::Upwind2}) {
        for (double a : {1.0, -1.0}) {
            INFO("scheme = " << static_cast<int>(scheme) << ", a = " << a);

            for (int j = j0; j <= j1; ++j) {
                for (int i = i0; i <= i1; ++i) {
                    INFO("i = " << i << ", j = " << j);
                    const auto [x, y] = position(g, s, i, j);

                    CHECK_THAT(derivative(fx, i, j, 1, 0, g.getDx(), a, scheme),
                               Catch::Matchers::WithinAbs(2.0 * x, 1e-10));
                    CHECK_THAT(derivative(fy, i, j, 0, 1, g.getDy(), a, scheme),
                               Catch::Matchers::WithinAbs(2.0 * y, 1e-10));
                }
            }
        }
    }
}

namespace {

// Manufactured velocities for the advection test (not divergence-free; we only test the operator)
double exactU(double x, double y) { return std::sin(pi * x) * std::cos(pi * y); }
double exactV(double x, double y) { return std::cos(pi * x) * std::sin(2.0 * pi * y); }

// Exact (u·∇)u and (u·∇)v for the fields above
double exactAdvU(double x, double y) {
    const double ux = pi * std::cos(pi * x) * std::cos(pi * y);
    const double uy = -pi * std::sin(pi * x) * std::sin(pi * y);
    return exactU(x, y) * ux + exactV(x, y) * uy;
}

double exactAdvV(double x, double y) {
    const double vx = -pi * std::sin(pi * x) * std::sin(2.0 * pi * y);
    const double vy = 2.0 * pi * std::cos(pi * x) * std::cos(2.0 * pi * y);
    return exactU(x, y) * vx + exactV(x, y) * vy;
}

// Max-norm error of advectionU (s = UFace) or advectionV (s = VFace) on an N x N grid with dx != dy
double maxAdvectionError(int N, AdvectionScheme scheme, Staggering s) {
    const Grid g(N, N, 2.0, 1.0);
    Field u(g);
    Field v(g);
    Field out(g);
    fillExact(u, g, Staggering::UFace, exactU);
    fillExact(v, g, Staggering::VFace, exactV);

    if (s == Staggering::UFace) {
        advectionU(u, v, g, scheme, out);
    } else {
        advectionV(u, v, g, scheme, out);
    }

    double err = 0.0;
    const auto [i0, i1, j0, j1] = interiorRange(g, s);
    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            const auto [x, y] = position(g, s, i, j);
            const double exact = (s == Staggering::UFace) ? exactAdvU(x, y) : exactAdvV(x, y);
            err = std::max(err, std::abs(out(i, j) - exact));
        }
    }
    return err;
}

} // namespace

TEST_CASE("Advection converges with second order", "[operators]") {
    for (AdvectionScheme scheme : {AdvectionScheme::Central, AdvectionScheme::Upwind2}) {
        for (Staggering s : {Staggering::UFace, Staggering::VFace}) {
            INFO("scheme = " << static_cast<int>(scheme) << ", staggering = " << static_cast<int>(s));

            const double e16 = maxAdvectionError(16, scheme, s);
            const double e32 = maxAdvectionError(32, scheme, s);
            const double e64 = maxAdvectionError(64, scheme, s);
            INFO("errors: " << e16 << " " << e32 << " " << e64);

            // One-sided: upwind converges slightly faster than 2 on coarse grids
            CHECK(std::log2(e16 / e32) > 1.85);
            CHECK(std::log2(e32 / e64) > 1.85);
        }
    }
}

TEST_CASE("cellCenteredVelocity is exact for linear fields", "[operators]") {
    const Grid g(5, 3, 2.0, 1.0);
    Field u(g);
    Field v(g);
    fillExact(u, g, Staggering::UFace, [](double x, double) { return 3.0 * x + 1.0; });
    fillExact(v, g, Staggering::VFace, [](double, double y) { return -2.0 * y; });

    Field uc(g, -999.0);
    Field vc(g, -999.0);
    cellCenteredVelocity(u, v, g, uc, vc);

    const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::Cell);
    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            INFO("i = " << i << ", j = " << j);
            const auto [x, y] = position(g, Staggering::Cell, i, j);
            CHECK_THAT(uc(i, j), Catch::Matchers::WithinAbs(3.0 * x + 1.0, 1e-12));
            CHECK_THAT(vc(i, j), Catch::Matchers::WithinAbs(-2.0 * y, 1e-12));
        }
    }
    CHECK(uc(0, 1) == -999.0);   // ghosts of the output are not written
}

TEST_CASE("vorticity is exact for linear velocity fields", "[operators]") {
    const Grid g(5, 3, 2.0, 1.0);   // dx != dy
    Field u(g);
    Field v(g);
    Field omega(g, -999.0);

    SECTION("solid-body rotation: w = 2") {
        fillExact(u, g, Staggering::UFace, [](double, double y) { return -(y - 0.5); });
        fillExact(v, g, Staggering::VFace, [](double x, double) { return x - 1.0; });
        vorticity(u, v, g, omega);
        for (int j = 1; j <= g.getNy(); ++j) {
            for (int i = 1; i <= g.getNx(); ++i) {
                INFO("i = " << i << ", j = " << j);
                CHECK_THAT(omega(i, j), Catch::Matchers::WithinAbs(2.0, 1e-12));
            }
        }
    }

    SECTION("simple shear u = y: w = -1") {
        fillExact(u, g, Staggering::UFace, [](double, double y) { return y; });
        fillExact(v, g, Staggering::VFace, [](double, double) { return 0.0; });
        vorticity(u, v, g, omega);
        for (int j = 1; j <= g.getNy(); ++j) {
            for (int i = 1; i <= g.getNx(); ++i) {
                INFO("i = " << i << ", j = " << j);
                CHECK_THAT(omega(i, j), Catch::Matchers::WithinAbs(-1.0, 1e-12));
            }
        }
    }

    SECTION("dv/dx and du/dy use the right spacing") {
        // u = 0, v = 4x  ->  w = 4 ; catches dx/dy swaps since dx = 0.4, dy = 1/3
        fillExact(u, g, Staggering::UFace, [](double, double) { return 0.0; });
        fillExact(v, g, Staggering::VFace, [](double x, double) { return 4.0 * x; });
        vorticity(u, v, g, omega);
        CHECK_THAT(omega(3, 2), Catch::Matchers::WithinAbs(4.0, 1e-12));
    }
}
