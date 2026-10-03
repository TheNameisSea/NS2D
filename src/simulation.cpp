#include "ns2d/simulation.h"

Simulation::Simulation(const Grid& grid_, const SimulationParams& params_) 
                        : grid(grid_), params(params_), u(grid), v(grid), p(grid),
                        uStar(grid), vStar(grid), work(grid){

    params.bc.applyNoSlipWall(u, v, grid);

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