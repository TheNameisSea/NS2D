#include "ns2d/io/summary.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

void writeSummary(const std::string& path, const RunSummary& summary) {
    nlohmann::json j;
    j["Re"] = summary.Re;
    j["nx"] = summary.nx;
    j["ny"] = summary.ny;
    j["converged"] = summary.converged;
    j["steps"] = summary.steps;
    j["time"] = summary.time;
    j["wallSeconds"] = summary.wallSeconds;
    j["finalChange"] = summary.finalChange;

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write summary file: " + path);
    }
    out << j.dump(2) << '\n';
    if (!out) {
        throw std::runtime_error("error while writing summary file: " + path);
    }
}
