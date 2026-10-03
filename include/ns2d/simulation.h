#pragma once
#include "ns2d/operators.h"
#include "ns2d/boundary.h"
#include "ns2d/projection.h"

struct SimulationParams {
    double Re{100.0};
    double dt{0.001};
    AdvectionScheme scheme{AdvectionScheme::Upwind2};
    BoundaryConditions bc{ .top = 1.0 };    // lid-driven cavity by default
};

class Simulation {
    private:
        Grid grid;                 // must be declared FIRST
        SimulationParams params;
        Field u, v, p;
        Field uStar, vStar;
        PredictorWork work;
    public:
        Simulation(const Grid& g, const SimulationParams& prm);

        void applyBoundaryConditions();
        void predict();

        const Field& getU() const;
        const Field& getV() const;
        const Field& getUStar() const;
        const Field& getVStar() const;

        Field& getU();

};

