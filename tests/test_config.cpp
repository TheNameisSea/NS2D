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
        "output":  { "printEvery": 7 }
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

    CHECK(cfg.steadyTol == 1e-5);
    CHECK(cfg.maxSteps == 1234);
    CHECK(cfg.printEvery == 7);
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

    // Case settings must get sensible (positive) defaults too
    CHECK(cfg.steadyTol > 0.0);
    CHECK(cfg.maxSteps > 0);
    CHECK(cfg.printEvery > 0);
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

    SECTION("the shipped cavity config loads") {
        const std::string path = std::string(NS2D_SOURCE_DIR) + "/configs/cavity_re100.json";
        const CaseConfig cfg = loadConfig(path);
        CHECK(cfg.grid.getNx() == 32);
        CHECK(cfg.grid.getNy() == 32);
        CHECK(cfg.params.Re == 100.0);
        CHECK(cfg.params.scheme == AdvectionScheme::Upwind2);
    }
}
