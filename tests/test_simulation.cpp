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

TEST_CASE("step() advances time by the adaptive dt", "[simulation]") {
    const Grid g(16, 16, 1.0, 1.0);
    const SimulationParams params{};            // adaptiveDt = true, cfl = 0.5, Re = 100
    Simulation sim(g, params);

    const double dtExpected = computeDt(sim.getU(), sim.getV(), g, params.Re, params.cfl);
    const double dt = sim.step().dt;

    CHECK(dt == dtExpected);
    CHECK(sim.getDt() == dt);
    CHECK(sim.getTime() == dt);

    const double dt2 = sim.step().dt;
    CHECK_THAT(sim.getTime(), Catch::Matchers::WithinRel(dt + dt2, 1e-14));
}

TEST_CASE("step() uses the fixed dt when adaptiveDt is off", "[simulation]") {
    const Grid g(8, 8, 1.0, 1.0);
    const SimulationParams params{.dt = 0.002, .adaptiveDt = false};
    Simulation sim(g, params);

    CHECK(sim.step().dt == 0.002);
    CHECK(sim.getDt() == 0.002);
    sim.step();
    CHECK_THAT(sim.getTime(), Catch::Matchers::WithinRel(0.004, 1e-14));
}

TEST_CASE("Cavity stays stable for 200 adaptive steps", "[simulation]") {
    const Grid g(16, 16, 1.0, 1.0);
    Simulation sim(g, SimulationParams{});

    for (int n = 0; n < 200; ++n) {
        sim.step();
    }

    const Field& u = sim.getU();
    const Field& v = sim.getV();
    const double uMax = maxAbs(u, interiorRange(g, Staggering::UFace));
    const double vMax = maxAbs(v, interiorRange(g, Staggering::VFace));
    INFO("t = " << sim.getTime() << ", max|u| = " << uMax << ", max|v| = " << vMax);

    CHECK(std::isfinite(uMax));
    CHECK(std::isfinite(vMax));
    CHECK(uMax < 1.5);                          // lid speed is 1
    CHECK(vMax < 1.5);
    CHECK(maxAbsDivergence(u, v, g) < 1e-6);   // still incompressible
}

TEST_CASE("step() reports Poisson statistics and time in StepInfo", "[simulation]") {
    const Grid g(16, 16, 1.0, 1.0);
    const SimulationParams params{};
    Simulation sim(g, params);

    const StepInfo first = sim.step();
    CHECK(first.time == first.dt);
    CHECK(first.poissonIterations > 0);
    CHECK(first.poissonIterations < params.poissonMaxIter);
    CHECK(first.poissonResidual < params.poissonTol);
    CHECK(first.change > 0.0);                  // the lid set the fluid in motion

    const StepInfo second = sim.step();
    CHECK_THAT(second.time, Catch::Matchers::WithinRel(first.dt + second.dt, 1e-14));
    CHECK(second.time == sim.getTime());
}

TEST_CASE("Cavity at Re = 100 reaches a steady state with a primary vortex", "[simulation]") {
    const Grid g(16, 16, 1.0, 1.0);
    Simulation sim(g, SimulationParams{});
    const double steadyTol = 1e-6;
    const int maxSteps = 5000;

    int steps = 0;
    StepInfo info{};
    for (steps = 1; steps <= maxSteps; ++steps) {
        info = sim.step();
        if (info.change < steadyTol) {
            break;
        }
    }
    INFO("steps = " << steps << ", t = " << info.time << ", change = " << info.change);

    REQUIRE(steps < maxSteps);                  // converged before the cap
    CHECK(info.change < steadyTol);

    // u along the vertical centerline x = 0.5 (u-faces at i = nx/2)
    const Field& u = sim.getU();
    const int ic = g.getNx() / 2;
    double uMin = 0.0;
    int jMin = 0;
    for (int j = 1; j <= g.getNy(); ++j) {
        if (u(ic, j) < uMin) {
            uMin = u(ic, j);
            jMin = j;
        }
    }
    INFO("min u on centerline = " << uMin << " at y = " << g.yp(jMin));

    // Primary vortex: flow returns (u < 0) in the lower half; Ghia et al. give min u ~ -0.21 at Re = 100
    CHECK(uMin < -0.1);
    CHECK(g.yp(jMin) < 0.5);
    CHECK(u(ic, g.getNy()) > 0.0);              // top row moves with the lid
}
