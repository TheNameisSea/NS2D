#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/boundary.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>
#include <algorithm>

TEST_CASE("bc test", "[bc]"){
    Grid test_grid(4, 3, 2.0, 1.0);
    Field u(test_grid, 1.5);
    Field v(test_grid, 1.5);

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(0.0, 10.0);

    std::generate(u.values().begin(), u.values().end(), [&]() { return dis(gen); });
    std::generate(v.values().begin(), v.values().end(), [&]() { return dis(gen); });

    Field u_copy{u};
    Field v_copy{v};

    BoundaryConditions bc{.top = 1.0};
    bc.applyNoSlipWall(u, v, test_grid);

    Field u_first_call{u};
    Field v_first_call{v};

    SECTION ("Apply boundary conditions") {
        // Normal velocities should be zero at the walls
        for (int j = 1; j <= test_grid.getNy(); ++j) {
            CHECK_THAT(u(0, j), Catch::Matchers::WithinAbs(0.0, 1e-14));
            CHECK_THAT(u(test_grid.getNx(), j), Catch::Matchers::WithinAbs(0.0, 1e-14));
        }
        for (int i = 1; i <= test_grid.getNx(); ++i) {
            CHECK_THAT(v(i, 0), Catch::Matchers::WithinAbs(0.0, 1e-14));
            CHECK_THAT(v(i, test_grid.getNy()), Catch::Matchers::WithinAbs(0.0, 1e-14));
        }

        // Tangential velocities should match the specified boundary conditions
        for (int i = 0; i <= test_grid.getNx(); ++i) {
            CHECK_THAT((u(i, test_grid.getNy()) + u(i, test_grid.getNy()+1)) / 2, Catch::Matchers::WithinAbs(1.0, 1e-14)); 
        }
        for (int i = 0; i <= test_grid.getNx(); ++i) {
            CHECK_THAT((u(i, 0) + u(i, 1)) / 2, Catch::Matchers::WithinAbs(0.0, 1e-14));
        }

        for (int j = 0; j <= test_grid.getNy(); ++j) {
            CHECK_THAT((v(test_grid.getNx(), j) + v(test_grid.getNx()+1, j)) / 2, Catch::Matchers::WithinAbs(0.0, 1e-14)); 
        }
        for (int j = 0; j <= test_grid.getNy(); ++j) {
            CHECK_THAT((v(0, j) + v(1, j)) / 2, Catch::Matchers::WithinAbs(0.0, 1e-14));
        }
    }

    SECTION ("Interior untouched, apply bc twice") {
        
        for (int j = 1; j <= test_grid.getNy(); ++j){
            for (int i = 1; i <= test_grid.getNx()-1; ++i){
                CHECK(u(i, j) == u_copy(i, j));

            }
        }

        for (int j = 1; j <= test_grid.getNy()-1; ++j){
            for (int i = 1; i <= test_grid.getNx(); ++i){
                CHECK(v(i, j) == v_copy(i, j));

            }
        }

        bc.applyNoSlipWall(u_first_call, v_first_call, test_grid);
        CHECK(u_first_call.values() == u.values());
        CHECK(v_first_call.values() == v.values());

    }

}
