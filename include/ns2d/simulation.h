#pragma once
#include "ns2d/operators.h"
#include "ns2d/boundary.h"
#include "ns2d/projection.h"
#include "ns2d/poisson.h"
#include <memory>

struct SimulationParams {
    double Re{100.0};
    double dt{0.001};
    AdvectionScheme scheme{AdvectionScheme::Upwind2};
    BoundaryConditions bc{ .top = 1.0 };    // lid-driven cavity by default
    int poissonMaxIter{100000};
    double poissonTol{1e-8};
};

class Simulation {
    private:
        Grid grid;                 // must be declared FIRST
        SimulationParams params;
        Field u, v, p;
        Field uStar, vStar;
        Field div, rhs;
        PredictorWork work;
        std::unique_ptr<PoissonSolver> solver;

    public:
        Simulation(const Grid& g, const SimulationParams& prm);
        
        SolveResult project();
        void applyBoundaryConditions();
        void predict();

        const Field& getU() const;
        const Field& getV() const;
        const Field& getUStar() const;
        const Field& getVStar() const;

        Field& getU();

};

