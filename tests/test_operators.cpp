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
