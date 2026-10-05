#include "ns2d/poisson.h"
#include <cmath>
#include <limits>
#include <utility>
#include <stdexcept>

namespace {
    void axpy(double a, const Field& x, Field& y, const Grid& g){
        
        const int nx = g.getNx();
        const int ny = g.getNy(); 

        for (int j = 1; j <= ny; ++j){
            for (int i = 1; i <= nx; ++i){
                y(i, j) += a * x(i, j);
            }
        }
    }
}

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

// Use plain sum
double interiorDot(const Field& a, const Field& b, const Grid& grid){
    
    const int nx = grid.getNx();
    const int ny = grid.getNy();

    double dotProduct{};

    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            dotProduct += a(i, j) * b(i, j);
        }
    }
    return dotProduct;

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

PCGSolver::PCGSolver(const Grid& g_, int maxIterations_)
                            : grid(g_), maxIterations(maxIterations_),
                            b(grid, 0.0), r(grid, 0.0), z(grid, 0.0), d(grid, 0.0), 
                            Ad(grid, 0.0), invDiag(grid, 0.0){

    const double dx = grid.getDx();
    const double dy = grid.getDy();
    const int nx = grid.getNx();
    const int ny = grid.getNy();  
    const double idx2 = 1.0/(dx*dx);  
    const double idy2 = 1.0/(dy*dy);

    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            const double diag = (i > 1  ? idx2 : 0) 
                        + (i < nx ? idx2 : 0) 
                        + (j > 1  ? idy2 : 0) 
                        + (j < ny ? idy2 : 0);
            invDiag(i,j) = 1.0 / diag;

        }
    }


}

SolveResult PCGSolver::solve(const Field& rhs, Field& p, double tol) {
    
    auto applyA = [&](Field& x, Field& out) {

        // Fill x's Neumann ghosts
        applyPressureNeumann(x, grid);
        // out = laplacian of x (Cell)
        laplacian(x, grid, Staggering::Cell, out);

        const int nx = grid.getNx();
        const int ny = grid.getNy(); 
        // Negate out over the interior
        for (int j = 1; j <= ny; ++j){
            for (int i = 1; i <= nx; ++i){
                out(i, j) = -out(i, j);
            }
        }
    };


    b = rhs;
    const int nx = grid.getNx();
    const int ny = grid.getNy(); 

    int iterations_done {maxIterations};

    const double bMean = interiorMean(b, grid);

    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            b(i,j) = -(b(i,j) - bMean);
        }
    }
    const double bNorm = interiorNorm(b, grid);

    if (bNorm == 0){
        p.fill(0.0);
        applyPressureNeumann(p, grid);
        return {0, 0};
    } 

    double res = std::numeric_limits<double>::infinity();

    applyA(p, Ad);
    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            r(i,j) = b(i,j) - Ad(i,j);
        }
    }

    // res = ‖r‖ / bNorm
    res = interiorNorm(r, grid) / bNorm;
    // if res < tol → skip the loop, iterations_done = 0
    if (res < tol){
        iterations_done = 0;
    }

    // z = invDiag · r     
    for (int j = 1; j <= ny; ++j){
        for (int i = 1; i <= nx; ++i){
            z(i,j) = invDiag(i,j) * r(i,j);
        }
    }

    d = z;

    double rz = interiorDot(r, z, grid);

    if (res >= tol) {
        for (int it = 1; it <= maxIterations; ++it){
            applyA(d, Ad);

            const double dAd = interiorDot(d, Ad, grid);
            if (dAd <= 0) {
                throw std::runtime_error("sign bug or missing ghosts");   // sign bug or missing ghosts
            }

            const double alpha = rz / dAd;
            // one interior loop:  p += alpha·d ;  r −= alpha·Ad
            axpy(alpha, d, p, grid);
            axpy(-alpha, Ad, r, grid);


            res = interiorNorm(r, grid) / bNorm;
            if (res < tol){
                iterations_done = it;
                break;
            }

            // one interior loop:  z = invDiag·r
            for (int j = 1; j <= ny; ++j){
                for (int i = 1; i <= nx; ++i){
                    z(i,j) = invDiag(i,j) * r(i,j);
                }
            }
            const double rzNew = interiorDot(r, z, grid);

            // one interior loop:  d = z + (rzNew/rz)·d
            for (int j = 1; j <= ny; ++j){
                for (int i = 1; i <= nx; ++i){
                    d(i,j) = z(i,j) + (rzNew/rz) * d(i,j);
                }
            }
            rz = rzNew;
        }
    }

    const double pMean = interiorMean(p, grid);
    for (int j = 1; j <= ny; ++j){
            for (int i = 1;  i<= nx; ++i){
                p(i, j) -= pMean;
            }
    }
    applyPressureNeumann(p, grid);
    return {iterations_done, res};

}

