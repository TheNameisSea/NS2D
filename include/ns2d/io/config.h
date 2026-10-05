#pragma once
#include "ns2d/grid.h"
#include "ns2d/simulation.h"
#include <nlohmann/json.hpp>
#include <string>

struct CaseConfig {
    Grid grid;
    SimulationParams params;
    double steadyTol{};
    int maxSteps{};
    int printEvery{};
    int vtkEvery{};                 // write a VTK file every N steps; 0 = only the final state
    std::string outputDir{"out"};   // folder for VTK output
};

CaseConfig parseConfig(const nlohmann::json& j);   // all the work; testable from a string
CaseConfig loadConfig(const std::string& path);    // open file → parse → parseConfig
