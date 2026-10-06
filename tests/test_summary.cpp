#include "ns2d/io/summary.h"
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

TEST_CASE("writeSummary writes every field as JSON", "[summary]") {
    const std::string path =
        (std::filesystem::temp_directory_path() / "ns2d_test_summary.json").string();
    const RunSummary written{
        .Re = 400.0,
        .nx = 64,
        .ny = 32,
        .converged = true,
        .steps = 19824,
        .time = 50.0883,
        .wallSeconds = 28.57,
        .finalChange = 9.9e-7,
    };

    writeSummary(path, written);

    std::ifstream in(path);
    REQUIRE(in);
    const nlohmann::json j = nlohmann::json::parse(in);

    CHECK(j.at("Re").get<double>() == written.Re);
    CHECK(j.at("nx").get<int>() == written.nx);
    CHECK(j.at("ny").get<int>() == written.ny);
    CHECK(j.at("converged").get<bool>() == written.converged);
    CHECK(j.at("steps").get<int>() == written.steps);
    CHECK(j.at("time").get<double>() == written.time);         // JSON keeps doubles exactly
    CHECK(j.at("wallSeconds").get<double>() == written.wallSeconds);
    CHECK(j.at("finalChange").get<double>() == written.finalChange);
    CHECK(j.size() == 8);                                         // nothing missing or extra

    in.close();
    std::filesystem::remove(path);
}

TEST_CASE("writeSummary records a run that did not converge", "[summary]") {
    const std::string path =
        (std::filesystem::temp_directory_path() / "ns2d_test_summary_nc.json").string();

    writeSummary(path, RunSummary{.Re = 1000.0, .nx = 8, .ny = 8, .converged = false, .steps = 10});

    std::ifstream in(path);
    REQUIRE(in);
    const nlohmann::json j = nlohmann::json::parse(in);
    CHECK(j.at("converged").get<bool>() == false);
    CHECK(j.at("steps").get<int>() == 10);

    in.close();
    std::filesystem::remove(path);
}

TEST_CASE("writeSummary throws when the file cannot be created", "[summary]") {
    CHECK_THROWS_AS(writeSummary("does/not/exist/summary.json", RunSummary{}), std::runtime_error);
}
