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
    PoissonSolverType poissonSolver{PoissonSolverType::PCG};
    int poissonMaxIter{100000};
    double poissonTol{1e-8};
    bool adaptiveDt{true};
    double cfl{0.25};
};


struct StepInfo {
    double dt{};
    double time{};
    int poissonIterations{};
    double poissonResidual{};
    double change{};
};

class Simulation {
    private:
        Grid grid;                 // must be declared FIRST
        SimulationParams params;
        Field u, v, p;
        Field uOld, vOld;
        Field uStar, vStar;
        Field div, rhs;
        PredictorWork work;
        std::unique_ptr<PoissonSolver> solver;
        double dt{};
        double time{0.0};

    public:
        Simulation(const Grid& g, const SimulationParams& prm);
        
        SolveResult project();
        void applyBoundaryConditions();
        void predict();


        const Field& getU() const;
        const Field& getV() const;
        const Field& getUStar() const;
        const Field& getVStar() const;
        const Field& getP() const;
        double getTime() const;
        double getDt() const;

        Field& getU();

        StepInfo step();

};




