#include "ns2d/io/config.h"
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

TEST_CASE("parseConfig reads every field of a full config", "[config]") {
    const auto j = nlohmann::json::parse(R"({
        "grid":    { "nx": 8, "ny": 4, "Lx": 2.0, "Ly": 1.0 },
        "Re": 400,
        "lidVelocity": 2.0,
        "scheme": "central",
        "cfl": 0.3,
        "poisson": { "solver": "jacobi", "tol": 1e-9, "maxIter": 5000 },
        "steady":  { "tol": 1e-5, "maxSteps": 1234 },
        "output":  { "printEvery": 7, "vtkEvery": 50, "directory": "results/run1" }
    })");

    const CaseConfig cfg = parseConfig(j);

    CHECK(cfg.grid.getNx() == 8);
    CHECK(cfg.grid.getNy() == 4);
    CHECK(cfg.grid.getLx() == 2.0);
    CHECK(cfg.grid.getLy() == 1.0);

    CHECK(cfg.params.Re == 400.0);
    CHECK(cfg.params.bc.top == 2.0);
    CHECK(cfg.params.bc.bottom == 0.0);
    CHECK(cfg.params.scheme == AdvectionScheme::Central);
    CHECK(cfg.params.cfl == 0.3);
    CHECK(cfg.params.poissonTol == 1e-9);
    CHECK(cfg.params.poissonMaxIter == 5000);
    CHECK(cfg.params.poissonSolver == PoissonSolverType::Jacobi);

    CHECK(cfg.steadyTol == 1e-5);
    CHECK(cfg.maxSteps == 1234);
    CHECK(cfg.printEvery == 7);
    CHECK(cfg.vtkEvery == 50);
    CHECK(cfg.outputDir == "results/run1");
}

TEST_CASE("parseConfig fills optional fields with defaults", "[config]") {
    // Only the required parts: grid and Re
    const auto j = nlohmann::json::parse(R"({
        "grid": { "nx": 16, "ny": 16, "Lx": 1.0, "Ly": 1.0 },
        "Re": 100
    })");

    const CaseConfig cfg = parseConfig(j);
    const SimulationParams defaults{};

    CHECK(cfg.params.Re == 100.0);
    CHECK(cfg.params.cfl == defaults.cfl);
    CHECK(cfg.params.scheme == defaults.scheme);
    CHECK(cfg.params.bc.top == defaults.bc.top);
    CHECK(cfg.params.poissonTol == defaults.poissonTol);
    CHECK(cfg.params.poissonMaxIter == defaults.poissonMaxIter);
    CHECK(cfg.params.poissonSolver == defaults.poissonSolver);
    CHECK(cfg.params.poissonSolver == PoissonSolverType::PCG);   // the fast one by default

    // Case settings must get sensible (positive) defaults too
    CHECK(cfg.steadyTol > 0.0);
    CHECK(cfg.maxSteps > 0);
    CHECK(cfg.printEvery > 0);
    CHECK(cfg.vtkEvery == 0);                   // default: only the final state
    CHECK(cfg.outputDir == "out");
}

TEST_CASE("parseConfig selects the Poisson solver", "[config]") {
    auto withPoisson = [](const std::string& poissonSection) {
        return nlohmann::json::parse(R"({
            "grid": { "nx": 16, "ny": 16, "Lx": 1.0, "Ly": 1.0 },
            "Re": 100,
            "poisson": )" + poissonSection + "}");
    };

    SECTION("\"pcg\"") {
        CHECK(parseConfig(withPoisson(R"({ "solver": "pcg" })")).params.poissonSolver
              == PoissonSolverType::PCG);
    }
    SECTION("\"jacobi\"") {
        CHECK(parseConfig(withPoisson(R"({ "solver": "jacobi" })")).params.poissonSolver
              == PoissonSolverType::Jacobi);
    }
    SECTION("poisson section without a solver key falls back to the default") {
        const CaseConfig cfg = parseConfig(withPoisson(R"({ "tol": 1e-9 })"));
        CHECK(cfg.params.poissonSolver == SimulationParams{}.poissonSolver);
        CHECK(cfg.params.poissonTol == 1e-9);
    }
}

TEST_CASE("parseConfig rejects bad input", "[config]") {
    SECTION("missing required Re") {
        const auto j = nlohmann::json::parse(R"({
            "grid": { "nx": 16, "ny": 16, "Lx": 1.0, "Ly": 1.0 }
        })");
        CHECK_THROWS(parseConfig(j));
    }

    SECTION("missing grid section") {
        const auto j = nlohmann::json::parse(R"({ "Re": 100 })");
        CHECK_THROWS(parseConfig(j));
    }

    SECTION("unknown advection scheme") {
        const auto j = nlohmann::json::parse(R"({
            "grid": { "nx": 16, "ny": 16, "Lx": 1.0, "Ly": 1.0 },
            "Re": 100,
            "scheme": "upwind3"
        })");
        CHECK_THROWS_AS(parseConfig(j), std::invalid_argument);
    }

    SECTION("unknown Poisson solver") {
        const auto j = nlohmann::json::parse(R"({
            "grid": { "nx": 16, "ny": 16, "Lx": 1.0, "Ly": 1.0 },
            "Re": 100,
            "poisson": { "solver": "magic" }
        })");
        CHECK_THROWS_AS(parseConfig(j), std::invalid_argument);
    }

    SECTION("wrong value type") {
        const auto j = nlohmann::json::parse(R"({
            "grid": { "nx": "sixteen", "ny": 16, "Lx": 1.0, "Ly": 1.0 },
            "Re": 100
        })");
        CHECK_THROWS(parseConfig(j));
    }

    SECTION("invalid grid size is caught by Grid") {
        const auto j = nlohmann::json::parse(R"({
            "grid": { "nx": 0, "ny": 16, "Lx": 1.0, "Ly": 1.0 },
            "Re": 100
        })");
        CHECK_THROWS_AS(parseConfig(j), std::invalid_argument);
    }
}

TEST_CASE("loadConfig reads files", "[config]") {
    SECTION("missing file throws runtime_error") {
        CHECK_THROWS_AS(loadConfig("does/not/exist.json"), std::runtime_error);
    }

    SECTION("the shipped cavity configs load with a fine enough grid") {
        // 1.15: one config per Reynolds number; the boundary layer thins with Re,
        // so the minimum grid grows (Re 1000 needs at least 128^2)
        struct Case { const char* file; double Re; int minN; };
        const Case cases[] = {
            {"cavity_re100.json", 100.0, 32},
            {"cavity_re400.json", 400.0, 64},
            {"cavity_re1000.json", 1000.0, 128},
        };
        for (const Case& c : cases) {
            INFO(c.file);
            const CaseConfig cfg = loadConfig(std::string(NS2D_SOURCE_DIR) + "/configs/" + c.file);
            CHECK(cfg.params.Re == c.Re);
            CHECK(cfg.grid.getNx() >= c.minN);
            CHECK(cfg.grid.getNy() == cfg.grid.getNx());
            CHECK(cfg.params.scheme == AdvectionScheme::Upwind2);
            CHECK(cfg.params.poissonSolver == PoissonSolverType::PCG);
            CHECK(cfg.params.cfl <= 0.25);      // upwind2 + forward Euler is unstable at 0.5 (64²)
            CHECK(cfg.outputDir.find(std::to_string(static_cast<int>(c.Re))) != std::string::npos);
        }
    }
}
