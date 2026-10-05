#pragma once
#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/operators.h"


struct SolveResult { int iterations{}; double residual{}; };

class PoissonSolver {
    public:
        virtual ~PoissonSolver() = default;
        virtual SolveResult solve(const Field& rhs, Field& p, double tol) = 0;
};

void applyPressureNeumann(Field& p, const Grid& grid);

double interiorMean(const Field& f,  const Grid& grid);

// Use RMS
double interiorNorm(const Field& f,  const Grid& grid);

// Use plain sum
double interiorDot(const Field& a, const Field& b, const Grid& grid);

class JacobiSolver : public PoissonSolver {
    public:
        JacobiSolver(const Grid& g, int maxIterations, int checkEvery = 10);
        SolveResult solve(const Field& rhs, Field& p, double tol) override;
    private:
        Grid grid;
        int maxIterations;
        int checkEvery;
        Field b;        // rhs with its mean removed
        Field pNew;     // second buffer for Jacobi
        Field r;        // Laplacian/residual work field
};

class PCGSolver : public PoissonSolver {
    public:
        PCGSolver(const Grid& g, int maxIterations);
        SolveResult solve(const Field& rhs, Field& p, double tol) override;
    private:
        Grid grid;
        int maxIterations;

        Field b;        // −(rhs − mean)
        Field r;        // Laplacian/residual work field
        Field z, d, Ad, invDiag;
};

