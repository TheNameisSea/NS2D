#include "ns2d/grid.h"
#include "ns2d/simulation.h"
#include <algorithm>
#include <cmath>
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

namespace {

double maxAbsDivergence(const Field& u, const Field& v, const Grid& g) {
    Field div(g, 0.0);
    divergence(u, v, g, div);
    double m = 0.0;
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            m = std::max(m, std::abs(div(i, j)));
        }
    }
    return m;
}

} // namespace

TEST_CASE("One full time step: project() makes the cavity flow divergence-free", "[simulation]") {
    const Grid g(16, 16, 1.0, 1.0);
    const SimulationParams params{};
    Simulation sim(g, params);

    sim.predict();
    const double divStar = maxAbsDivergence(sim.getUStar(), sim.getVStar(), g);

    const SolveResult res = sim.project();
    const double divNew = maxAbsDivergence(sim.getU(), sim.getV(), g);
    INFO("Poisson: iterations = " << res.iterations << ", residual = " << res.residual);
    INFO("max|div u*| = " << divStar << ", max|div u| = " << divNew);

    CHECK(divStar > 1e-3);                         // the predictor really made a divergent field (~0.08 here)
    CHECK(res.residual < params.poissonTol);       // the solver converged
    CHECK(res.iterations < params.poissonMaxIter);
    CHECK(divNew / divStar < 1e-6);

    SECTION("the lid BC still holds after projection") {
        const Field& u = sim.getU();
        const int ny = g.getNy();
        for (int i = 1; i <= g.getNx() - 1; ++i) {
            INFO("i = " << i);
            CHECK_THAT(0.5 * (u(i, ny) + u(i, ny + 1)), Catch::Matchers::WithinAbs(1.0, 1e-12));
        }
    }

    SECTION("the flow has started moving under the lid") {
        // After one step the top interior row is pulled in the lid direction
        const Field& u = sim.getU();
        double topRowSum = 0.0;
        for (int i = 1; i <= g.getNx() - 1; ++i) {
            topRowSum += u(i, g.getNy());
        }
        CHECK(topRowSum > 0.0);
    }
}
