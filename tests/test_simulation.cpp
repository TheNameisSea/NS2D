#include "ns2d/grid.h"
#include "ns2d/simulation.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Simulation applies the lid BC on construction", "[simulation]") {
    const Grid g(8, 8, 1.0, 1.0);
    const Simulation sim(g, SimulationParams{});   // default: lid u = 1, fluid at rest
    const Field& u = sim.getU();

    // Lid ghost: average of u(i, ny) = 0 and u(i, ny+1) must be 1
    for (int i = 0; i <= g.getNx(); ++i) {
        INFO("i = " << i);
        CHECK_THAT(u(i, g.getNy() + 1), Catch::Matchers::WithinAbs(2.0, 1e-14));
    }
}

TEST_CASE("One predictor step: the lid drags the top row, nothing else moves", "[simulation]") {
    const Grid g(8, 8, 1.0, 1.0);
    Simulation sim(g, SimulationParams{});
    sim.predict();

    const Field& uStar = sim.getUStar();
    const int nx = g.getNx();
    const int ny = g.getNy();

    for (int i = 1; i <= nx - 1; ++i) {
        INFO("i = " << i);
        // Diffusion from the moving lid pulls the top interior row forward
        CHECK(uStar(i, ny) > 0.0);
        // In one explicit step the 5-point stencil only reaches one row
        CHECK(uStar(i, ny - 1) == 0.0);
    }

    // Walls stay at rest
    for (int j = 1; j <= ny; ++j) {
        INFO("j = " << j);
        CHECK(uStar(0, j) == 0.0);
        CHECK(uStar(nx, j) == 0.0);
    }
}
