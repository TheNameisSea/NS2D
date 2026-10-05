#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/operators.h"
#include "ns2d/poisson.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <algorithm>
#include <random>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>

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

// Manufactured solution with dp/dn = 0 on all walls of the unit square.
// Note: cos(πx)cos(πy) is an eigenvector of the discrete Neumann Laplacian, so plain CG
// would solve it in one iteration. PCG tests therefore also use a random right-hand side.
double exactP(double x, double y) {
    const double pi = 3.14159265358979323846;
    return std::cos(pi * x) * std::cos(pi * y);
}

double exactRhs(double x, double y) {
    const double pi = 3.14159265358979323846;
    return -2.0 * pi * pi * exactP(x, y);
}

enum class SolverKind { Jacobi, PCG };

const char* name(SolverKind kind) {
    return kind == SolverKind::Jacobi ? "Jacobi" : "PCG";
}

std::unique_ptr<PoissonSolver> makeSolver(SolverKind kind, const Grid& g, int maxIter) {
    switch (kind) {
        case SolverKind::Jacobi: return std::make_unique<JacobiSolver>(g, maxIter);
        case SolverKind::PCG:    return std::make_unique<PCGSolver>(g, maxIter);
    }
    return nullptr;
}

struct PoissonRun {
    double maxError{};
    double pMean{};
    SolveResult result{};
};

// Solve the manufactured problem on an N x N unit square through the PoissonSolver interface
PoissonRun runManufactured(SolverKind kind, int N, double tol) {
    const Grid g(N, N, 1.0, 1.0);
    Field rhs(g);
    Field p(g, 0.0);
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            rhs(i, j) = exactRhs(g.xp(i), g.yp(j));
        }
    }

    const auto solver = makeSolver(kind, g, 200000);
    const SolveResult res = solver->solve(rhs, p, tol);

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

// Random right-hand side in [-1, 1] (not mean-free: the solver must handle that)
Field randomRhs(const Grid& g, unsigned seed) {
    Field rhs(g);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            rhs(i, j) = dis(gen);
        }
    }
    return rhs;
}

// Residual computed from scratch: ||∇²p − (rhs − mean)|| / ||rhs − mean|| (RMS norms)
double trueResidual(const Field& rhs, Field p, const Grid& g) {
    applyPressureNeumann(p, g);
    Field lap(g);
    laplacian(p, g, Staggering::Cell, lap);

    const double mean = interiorMean(rhs, g);
    Field b(g);
    Field r(g);
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            b(i, j) = rhs(i, j) - mean;
            r(i, j) = lap(i, j) - b(i, j);
        }
    }
    return interiorNorm(r, g) / interiorNorm(b, g);
}

} // namespace

// ---------------------------------------------------------------------------
// Tests shared by every solver (GENERATE runs the test case once per solver)
// ---------------------------------------------------------------------------

TEST_CASE("Poisson solvers converge to the manufactured solution", "[poisson]") {
    const SolverKind kind = GENERATE(SolverKind::Jacobi, SolverKind::PCG);
    const double tol = 1e-10;
    const PoissonRun run = runManufactured(kind, 16, tol);
    INFO(name(kind) << ": iterations = " << run.result.iterations
                    << ", residual = " << run.result.residual);

    CHECK(run.result.residual < tol);
    CHECK(run.result.iterations > 0);
    CHECK(run.result.iterations < 200000);              // converged, did not hit the cap
    CHECK_THAT(run.pMean, Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK(run.maxError < 1e-2);                         // discretization error at N = 16
}

TEST_CASE("Poisson solutions converge with second order in space", "[poisson]") {
    const SolverKind kind = GENERATE(SolverKind::Jacobi, SolverKind::PCG);
    const double tol = 1e-10;
    const double e8 = runManufactured(kind, 8, tol).maxError;
    const double e16 = runManufactured(kind, 16, tol).maxError;
    const double e32 = runManufactured(kind, 32, tol).maxError;
    INFO(name(kind) << " errors: " << e8 << " " << e16 << " " << e32);

    CHECK(std::log2(e8 / e16) > 1.9);
    CHECK(std::log2(e16 / e32) > 1.9);
}

TEST_CASE("Poisson solvers with zero right-hand side return p = 0 immediately", "[poisson]") {
    const SolverKind kind = GENERATE(SolverKind::Jacobi, SolverKind::PCG);
    INFO(name(kind));
    const Grid g(6, 4, 2.0, 1.0);
    const Field rhs(g, 0.0);
    Field p(g, 3.0);
    const auto solver = makeSolver(kind, g, 1000);

    const SolveResult res = solver->solve(rhs, p, 1e-10);

    CHECK(res.iterations == 0);
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            CHECK(p(i, j) == 0.0);
        }
    }
}

TEST_CASE("Poisson solvers remove a constant offset in the right-hand side", "[poisson]") {
    // rhs + c is not solvable with pure Neumann BCs; the solver must subtract the mean
    const SolverKind kind = GENERATE(SolverKind::Jacobi, SolverKind::PCG);
    const int N = 16;
    const Grid g(N, N, 1.0, 1.0);
    Field rhs(g);
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            rhs(i, j) = exactRhs(g.xp(i), g.yp(j)) + 0.5;
        }
    }
    Field p(g, 0.0);
    const auto solver = makeSolver(kind, g, 200000);
    const SolveResult res = solver->solve(rhs, p, 1e-10);
    INFO(name(kind) << ": iterations = " << res.iterations << ", residual = " << res.residual);

    CHECK(res.residual < 1e-10);
    CHECK(res.iterations < 200000);
}

TEST_CASE("Poisson solvers solve a random right-hand side on a non-square grid", "[poisson]") {
    // A random rhs mixes all modes, so it is a much harder test than one eigenvector
    const SolverKind kind = GENERATE(SolverKind::Jacobi, SolverKind::PCG);
    const Grid g(12, 8, 1.5, 1.0);
    const Field rhs = randomRhs(g, 7);
    Field p(g, 0.0);
    const auto solver = makeSolver(kind, g, 200000);

    const SolveResult res = solver->solve(rhs, p, 1e-10);
    const double trueRes = trueResidual(rhs, p, g);
    INFO(name(kind) << ": iterations = " << res.iterations << ", reported = " << res.residual
                    << ", true = " << trueRes);

    CHECK(res.residual < 1e-10);
    CHECK(trueRes < 1e-9);                              // the answer really solves the system
    CHECK_THAT(interiorMean(p, g), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

// ---------------------------------------------------------------------------
// PCG-specific tests
// ---------------------------------------------------------------------------

TEST_CASE("interiorDot sums products over interior cells only", "[poisson]") {
    const Grid g(3, 2, 1.0, 1.0);
    Field a(g, 100.0);                                  // ghosts stay 100: must be ignored
    Field b(g, 100.0);
    for (int j = 1; j <= 2; ++j) {
        for (int i = 1; i <= 3; ++i) {
            a(i, j) = i;
            b(i, j) = j;
        }
    }
    // sum_{i=1..3} sum_{j=1..2} i * j = 6 * 3 = 18 (plain sum, not RMS)
    CHECK(interiorDot(a, b, g) == 18.0);
    CHECK(interiorDot(a, b, g) == interiorDot(b, a, g));
}

TEST_CASE("PCG needs far fewer iterations than Jacobi", "[poisson]") {
    const int jacobiIt = runManufactured(SolverKind::Jacobi, 16, 1e-10).result.iterations;
    const int pcgIt = runManufactured(SolverKind::PCG, 16, 1e-10).result.iterations;
    INFO("Jacobi " << jacobiIt << ", PCG " << pcgIt);   // reference: 1190 vs 30

    CHECK(pcgIt <= 60);
    CHECK(10 * pcgIt < jacobiIt);
}

TEST_CASE("PCG iterations grow about linearly with N", "[poisson]") {
    // Random rhs, unit square. Reference: N = 16 / 32 → 68 / 141 iterations
    auto iterations = [](int N) {
        const Grid g(N, N, 1.0, 1.0);
        const Field rhs = randomRhs(g, 1);
        Field p(g, 0.0);
        PCGSolver pcg(g, 100000);
        const SolveResult res = pcg.solve(rhs, p, 1e-10);
        CHECK(res.residual < 1e-10);
        return res.iterations;
    };
    const int it16 = iterations(16);
    const int it32 = iterations(32);
    INFO("N=16: " << it16 << ", N=32: " << it32);

    CHECK(it16 <= 120);
    CHECK(it32 < 3 * it16);                             // ~2x for O(N); Jacobi would be ~4x
}

TEST_CASE("PCG reports the true residual", "[poisson]") {
    const Grid g(16, 16, 1.0, 1.0);
    const Field rhs = randomRhs(g, 3);
    Field p(g, 0.0);
    PCGSolver pcg(g, 100000);

    const SolveResult res = pcg.solve(rhs, p, 1e-10);
    const double trueRes = trueResidual(rhs, p, g);
    INFO("reported " << res.residual << ", true " << trueRes);

    // The recursive residual r -= alpha*Ad must stay in step with b − A p
    CHECK_THAT(res.residual, Catch::Matchers::WithinRel(trueRes, 1e-3));
}

TEST_CASE("PCG warm start from a converged solution needs no iterations", "[poisson]") {
    const Grid g(16, 16, 1.0, 1.0);
    const Field rhs = randomRhs(g, 5);
    Field p(g, 0.0);
    PCGSolver pcg(g, 100000);

    const SolveResult first = pcg.solve(rhs, p, 1e-10);
    REQUIRE(first.residual < 1e-10);
    const Field converged{p};

    const SolveResult second = pcg.solve(rhs, p, 1e-10);

    CHECK(second.iterations == 0);
    CHECK(second.residual < 1e-10);
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            CHECK_THAT(p(i, j), Catch::Matchers::WithinAbs(converged(i, j), 1e-12));
        }
    }
}

TEST_CASE("PCG stops at the iteration cap and reports non-convergence", "[poisson]") {
    const Grid g(16, 16, 1.0, 1.0);
    const Field rhs = randomRhs(g, 9);
    Field p(g, 0.0);
    PCGSolver pcg(g, 5);

    const SolveResult res = pcg.solve(rhs, p, 1e-10);

    CHECK(res.iterations == 5);
    CHECK(res.residual > 1e-10);
    CHECK(std::isfinite(res.residual));
}
