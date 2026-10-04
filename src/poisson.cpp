#include "ns2d/poisson.h"
#include <cmath>
#include <limits>
#include <utility>

void applyPressureNeumann(Field& p, const Grid& grid){

    const int nx = grid.getNx();
    const int ny = grid.getNy();

    for (int j = 1; j <= ny; ++j){
        p(0, j) = p(1, j);
        p(nx+1, j) = p(nx, j);
            
    }

    for (int i = 1; i <= nx; ++i){
        p(i, 0) = p(i, 1);
        p(i, ny+1) = p(i, ny);
    }

}

double interiorMean(const Field& f,  const Grid& grid){

    const int nx = grid.getNx();
    const int ny = grid.getNy();

    double sum{};

    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            sum += f(i, j);
        }
    }
    return sum / (nx * ny);

}

// Use RMS
double interiorNorm(const Field& f, const Grid& grid){

    const int nx = grid.getNx();
    const int ny = grid.getNy();

    double sum_of_square{};

    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            sum_of_square += f(i, j) * f(i, j);
        }
    }
    return std::sqrt(sum_of_square / (nx * ny));

}   


JacobiSolver::JacobiSolver(const Grid& g_, int maxIterations_, int checkEvery_)
                            : grid(g_), maxIterations(maxIterations_), checkEvery(checkEvery_),
                            b(grid, 0.0), pNew(grid, 0.0), r(grid, 0.0){

}

SolveResult JacobiSolver::solve(const Field& rhs, Field& p, double tol) {
    
    b = rhs;
    const int nx = grid.getNx();
    const int ny = grid.getNy(); 
    const double dx = grid.getDx();
    const double dy = grid.getDy();
    int iterations_done {maxIterations};

    const double bMean = interiorMean(b, grid);

    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            b(i, j) -= bMean;
        }
    }
    const double bNorm = interiorNorm(b, grid);

    if (bNorm == 0){
        p.fill(0.0);
        applyPressureNeumann(p, grid);
        return {0, 0};
    } 

    const double denom = 2/(dx*dx) + 2/(dy*dy);  
    const double idx2 = 1/(dx*dx);  
    const double idy2 = 1/(dy*dy);

    double res = std::numeric_limits<double>::infinity();

    const auto [i0, i1, j0, j1] = interiorRange(grid, Staggering::Cell);

    for (int it = 1; it <= maxIterations; ++it){
        applyPressureNeumann(p, grid);
        for (int j = j0; j <= j1; ++j){
            for (int i = i0;  i<= i1; ++i){
                pNew(i,j) = ( (p(i+1,j)+ p(i-1,j)) * idx2 + (p(i,j+1) + p(i,j-1)) * idy2 - b(i,j) ) / denom;
            }
        }

        std::swap(p, pNew);
        if (it % checkEvery == 0){
            applyPressureNeumann(p, grid);
            laplacian(p, grid, Staggering::Cell, r);
            applyPressureNeumann(p, grid);
            for (int j = j0; j <= j1; ++j){
                for (int i = i0;  i<= i1; ++i){
                    r(i,j) = b(i,j) - r(i,j);
                }
            }
            res = interiorNorm(r, grid) / bNorm;
            if (res < tol){
                iterations_done = it;
                break;
            }  
        }
    }

    const double pMean = interiorMean(p, grid);
    for (int j = j0; j <= j1; ++j){
            for (int i = i0;  i<= i1; ++i){
                p(i, j) -= pMean;
            }
    }
    applyPressureNeumann(p, grid);
    return {iterations_done, res};

}
