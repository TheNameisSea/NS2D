// Lid-driven cavity: run to steady state from a JSON config.
// Usage: cavity <config.json>
// Exit codes: 0 = steady state reached, 1 = bad input / error, 2 = not converged.

#include "ns2d/io/config.h"
#include "ns2d/simulation.h"
#include <chrono>
#include <exception>
#include <format>
#include <iostream>

namespace {

void printHeader() {
    std::cout << std::format("{:>8} {:>10} {:>10} {:>8} {:>10} {:>10}\n",
                             "step", "t", "dt", "poisson", "residual", "change");
}

void printRow(int step, const StepInfo& info) {
    std::cout << std::format("{:>8} {:>10.4f} {:>10.3e} {:>8} {:>10.3e} {:>10.3e}\n",
                             step, info.time, info.dt, info.poissonIterations,
                             info.poissonResidual, info.change);
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: cavity <config.json>\n";
        return 1;
    }

    try {
        const CaseConfig cfg = loadConfig(argv[1]);
        Simulation sim(cfg.grid, cfg.params);

        std::cout << std::format("Cavity {}x{}, Re = {}, steady tol = {:.1e}, max steps = {}\n",
                                 cfg.grid.getNx(), cfg.grid.getNy(), cfg.params.Re,
                                 cfg.steadyTol, cfg.maxSteps);
        printHeader();

        const auto start = std::chrono::steady_clock::now();
        bool converged = false;
        int steps = 0;
        StepInfo info{};

        for (steps = 1; steps <= cfg.maxSteps; ++steps) {
            info = sim.step();
            converged = info.change < cfg.steadyTol;
            if (steps % cfg.printEvery == 0 || converged) {
                printRow(steps, info);
            }
            if (converged) {
                break;
            }
        }

        const double wallSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        if (!converged) {
            std::cerr << std::format("Not converged after {} steps (change = {:.3e} > {:.1e}), "
                                     "wall time {:.2f} s\n",
                                     cfg.maxSteps, info.change, cfg.steadyTol, wallSeconds);
            return 2;
        }

        std::cout << std::format("Steady state at step {}, t = {:.4f}, wall time {:.2f} s\n",
                                 steps, info.time, wallSeconds);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
