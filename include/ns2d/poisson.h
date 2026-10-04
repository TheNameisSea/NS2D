#pragma once
#include "ns2d/grid.h"
#include "ns2d/field.h"


struct SolveResult { int iterations{}; double residual{}; };

class PoissonSolver {
    public:
        virtual ~PoissonSolver() = default;
        virtual SolveResult solve(const Field& rhs, Field& p, double tol) = 0;
};

void applyPressureNeumann(Field& p, const Grid& grid);
