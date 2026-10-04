#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/poisson.h"
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <random>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

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

namespace {

// Manufactured solution with dp/dn = 0 on all walls of the unit square
double exactP(double x, double y) {
    const double pi = 3.14159265358979323846;
    return std::cos(pi * x) * std::cos(pi * y);
}

double exactRhs(double x, double y) {
    const double pi = 3.14159265358979323846;
    return -2.0 * pi * pi * exactP(x, y);
}

struct PoissonRun {
    double maxError{};
    double pMean{};
    SolveResult result{};
};

// Solve on an N x N unit square through the PoissonSolver interface
PoissonRun runJacobi(int N, double tol) {
    const Grid g(N, N, 1.0, 1.0);
    Field rhs(g);
    Field p(g, 0.0);
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            rhs(i, j) = exactRhs(g.xp(i), g.yp(j));
        }
    }

    JacobiSolver jacobi(g, 200000);
    PoissonSolver& solver = jacobi;   // use it only through the interface
    const SolveResult res = solver.solve(rhs, p, tol);

    // The solver fixes mean(p) = 0, so compare against the exact solution with its mean removed
    double exactMean = 0.0;
    double pSum = 0.0;
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            exactMean += exactP(g.xp(i), g.yp(j));
            pSum += p(i, j);
        }
    }
    exactMean /= N * N;

    double err = 0.0;
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            const double exact = exactP(g.xp(i), g.yp(j)) - exactMean;
            err = std::max(err, std::abs(p(i, j) - exact));
        }
    }
    return PoissonRun{err, pSum / (N * N), res};
}

} // namespace

TEST_CASE("Jacobi converges to the manufactured solution", "[poisson]") {
    const double tol = 1e-10;
    const PoissonRun run = runJacobi(16, tol);
    INFO("iterations = " << run.result.iterations << ", residual = " << run.result.residual);

    CHECK(run.result.residual < tol);
    CHECK(run.result.iterations > 0);
    CHECK(run.result.iterations < 200000);              // converged, did not hit the cap
    CHECK_THAT(run.pMean, Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK(run.maxError < 1e-2);                         // discretization error at N = 16
}

TEST_CASE("Jacobi solution converges with second order in space", "[poisson]") {
    const double tol = 1e-10;
    const double e8 = runJacobi(8, tol).maxError;
    const double e16 = runJacobi(16, tol).maxError;
    const double e32 = runJacobi(32, tol).maxError;
    INFO("errors: " << e8 << " " << e16 << " " << e32);

    CHECK(std::log2(e8 / e16) > 1.9);
    CHECK(std::log2(e16 / e32) > 1.9);
}

TEST_CASE("Jacobi with zero right-hand side returns p = 0 immediately", "[poisson]") {
    const Grid g(6, 4, 2.0, 1.0);
    const Field rhs(g, 0.0);
    Field p(g, 3.0);
    JacobiSolver jacobi(g, 1000);

    const SolveResult res = jacobi.solve(rhs, p, 1e-10);

    CHECK(res.iterations == 0);
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            CHECK(p(i, j) == 0.0);
        }
    }
}

TEST_CASE("Jacobi removes a constant offset in the right-hand side", "[poisson]") {
    // rhs + c is not solvable with pure Neumann BCs; the solver must subtract the mean
    const int N = 16;
    const Grid g(N, N, 1.0, 1.0);
    Field rhs(g);
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            rhs(i, j) = exactRhs(g.xp(i), g.yp(j)) + 0.5;
        }
    }
    Field p(g, 0.0);
    JacobiSolver jacobi(g, 200000);
    const SolveResult res = jacobi.solve(rhs, p, 1e-10);
    INFO("iterations = " << res.iterations << ", residual = " << res.residual);

    CHECK(res.residual < 1e-10);
    CHECK(res.iterations < 200000);
}
