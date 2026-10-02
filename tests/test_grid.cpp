#include "ns2d/grid.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <stdexcept>

TEST_CASE("grid test", "[grid]"){
    Grid test_grid(4, 3, 2.0, 1.0);
    CHECK(test_grid.getDx() == Catch::Approx(0.5));
    CHECK(test_grid.getDy() == Catch::Approx(1.0 / 3.0));
    CHECK(test_grid.getNi() == 6);
    CHECK(test_grid.getNj() == 5);
    CHECK(test_grid.getLx() == Catch::Approx(2.0));
    CHECK(test_grid.getLy() == Catch::Approx(1.0));
    CHECK(test_grid.size() == 30);
    CHECK(test_grid.xp(1) == Catch::Approx(0.25));
    CHECK(test_grid.yp(1) == Catch::Approx(1.0 / 6.0));
    CHECK(test_grid.xu(1) == Catch::Approx(0.5));
    CHECK(test_grid.yv(1) == Catch::Approx(1.0 / 3.0));
    CHECK_THROWS_AS(Grid(0, 3, 1.0, 1.0), std::invalid_argument);

    CHECK(test_grid.xu(0) == Catch::Approx(0.0));
    CHECK(test_grid.xu(test_grid.getNx()) == Catch::Approx(test_grid.getLx()));
    CHECK(test_grid.yv(test_grid.getNy()) == Catch::Approx(test_grid.getLy()));

    CHECK(test_grid.xp(0) == Catch::Approx(test_grid.getDx() * -0.5));
    CHECK_THROWS_AS(Grid(4, 3, -1.0, 1.0), std::invalid_argument);

}