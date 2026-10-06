#pragma once
#include <string>

// End-of-run facts for one case, written next to the VTK output so that scripts
// (run_cavity.sh, the 1.16 plots, the Phase-5 sweep) don't have to parse console output.
struct RunSummary {
    double Re{};
    int nx{};
    int ny{};
    bool converged{};
    int steps{};            // steps actually taken
    double time{};          // simulated time at the end
    double wallSeconds{};
    double finalChange{};   // max|du|/dt of the last step
};

// Write the summary as pretty-printed JSON. Throws std::runtime_error if the file cannot be written.
void writeSummary(const std::string& path, const RunSummary& summary);
