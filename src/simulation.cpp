#include "ns2d/simulation.h"

Simulation::Simulation(const Grid& grid_, const SimulationParams& params_) 
                        : grid(grid_), params(params_), u(grid), v(grid), p(grid),
                        uStar(grid), vStar(grid), div(grid), rhs(grid), work(grid), 
                        solver(std::make_unique<JacobiSolver>(grid, params.poissonMaxIter)){

    applyBoundaryConditions();

}

void Simulation::applyBoundaryConditions(){
    params.bc.applyNoSlipWall(u, v, grid);
}

void Simulation::predict() {
    applyBoundaryConditions();
    predictor(u, v, grid, params.dt, params.Re, params.scheme, work, uStar, vStar);
}

const Field& Simulation::getU() const{
    return u;
}
const Field& Simulation::getV() const{
    return v;
}
const Field& Simulation::getUStar() const{
    return uStar;
}
const Field& Simulation::getVStar() const{
    return vStar;
}

Field& Simulation::getU() {
    return u;
}

SolveResult Simulation::project(){

    const double one_over_dt{1.0/params.dt};

    // 1. Divergence of the predicted velocity, at cell centers
    divergence(uStar, vStar, grid, div);

    // 2. Right-hand side of the pressure equation
    const auto[i0, i1, j0, j1] = interiorRange(grid, Staggering::Cell);

    for (int j = j0; j <= j1; ++j){
        for (int i = i0; i <= i1; ++i){
            rhs(i, j) = div(i, j) * one_over_dt;
        }
    }

    // 3. Solve ∇²p = rhs, through the base-class pointer
    const SolveResult res = solver->solve(rhs, p, params.poissonTol);

    // 4. Correct the velocity: u = u* − dt ∇p
    corrector(uStar, vStar, p, grid, params.dt, u, v);

    // 5. Make the new velocity satisfy the BCs (tangential ghosts, lid)
    applyBoundaryConditions();

    return res;

}

