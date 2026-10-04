#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/poisson.h"
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <random>

TEST_CASE("applyPressureNeumann copies interior neighbours into ghosts", "[poisson]") {
    const Grid g(5, 3, 2.0, 1.0);
    const int nx = g.getNx();
    const int ny = g.getNy();

    Field p(g);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-10.0, 10.0);
    std::generate(p.values().begin(), p.values().end(), [&]() { return dis(gen); });
    const Field before{p};

    applyPressureNeumann(p, g);

    SECTION("ghosts equal their interior neighbour (dp/dn = 0)") {
        for (int j = 1; j <= ny; ++j) {
            INFO("j = " << j);
            CHECK(p(0, j) == p(1, j));            // left
            CHECK(p(nx + 1, j) == p(nx, j));      // right
        }
        for (int i = 1; i <= nx; ++i) {
            INFO("i = " << i);
            CHECK(p(i, 0) == p(i, 1));            // bottom
            CHECK(p(i, ny + 1) == p(i, ny));      // top
        }
    }

    SECTION("interior is untouched") {
        for (int j = 1; j <= ny; ++j) {
            for (int i = 1; i <= nx; ++i) {
                INFO("i = " << i << ", j = " << j);
                CHECK(p(i, j) == before(i, j));
            }
        }
    }

    SECTION("applying twice changes nothing") {
        const Field once{p};
        applyPressureNeumann(p, g);
        CHECK(p.values() == once.values());
    }
}
