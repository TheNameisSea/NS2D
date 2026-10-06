#include "ns2d/simulation.h"
#include <algorithm>
#include <cmath>

Simulation::Simulation(const Grid& grid_, const SimulationParams& params_) 
                        : grid(grid_), params(params_), u(grid), v(grid), p(grid),
                        uOld(grid), vOld(grid), uStar(grid), vStar(grid), div(grid), rhs(grid), work(grid), 
                        solver(makePoissonSolver(params.poissonSolver, grid, params.poissonMaxIter)), dt(params.dt){

    applyBoundaryConditions();

}

void Simulation::applyBoundaryConditions(){
    params.bc.applyNoSlipWall(u, v, grid);
}

void Simulation::predict() {
    applyBoundaryConditions();
    predictor(u, v, grid, dt, params.Re, params.scheme, work, uStar, vStar);
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
const Field& Simulation::getP() const{
    return p;
}

double Simulation::getTime() const {
    return time;
}    

double Simulation::getDt() const {
    return dt;
}

Field& Simulation::getU() {
    return u;
}

SolveResult Simulation::project(){

    const double one_over_dt{1.0/dt};

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
    corrector(uStar, vStar, p, grid, dt, u, v);

    // 5. Make the new velocity satisfy the BCs (tangential ghosts, lid)
    applyBoundaryConditions();

    return res;

}

StepInfo Simulation::step(){

    uOld = u; 
    vOld = v;

    dt = (params.adaptiveDt) ? computeDt(u, v, grid, params.Re, params.cfl) : params.dt;
    predict();
    const SolveResult res = project();
    time += dt;
    
    double max_change_u{0.0};
    {
        const auto[i0, i1, j0, j1] = interiorRange(grid, Staggering::UFace);

        for (int j = j0; j <= j1; ++j){
            for (int i = i0; i <= i1; ++i){
                max_change_u = std::max(max_change_u, std::abs(u(i,j) - uOld(i,j))/dt);
            }
        }   
    } 
    double max_change_v{0.0};
    {
        const auto[i0, i1, j0, j1] = interiorRange(grid, Staggering::VFace);

        for (int j = j0; j <= j1; ++j){
            for (int i = i0; i <= i1; ++i){
                max_change_v = std::max(max_change_v, std::abs(v(i,j) - vOld(i,j))/dt);
            }
        }   
    } 
    double change = std::max(max_change_u, max_change_v);

    return StepInfo{.dt = dt, .time=time, .poissonIterations=res.iterations, .poissonResidual=res.residual, .change=change};


}



