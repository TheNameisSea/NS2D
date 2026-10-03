#include "ns2d/operators.h"
#include "ns2d/grid.h"
#include <cassert>

IndexRange interiorRange(const Grid& g, Staggering s){

    const int nx = g.getNx();
    const int ny = g.getNy();

    switch (s)
    {
    case Staggering::Cell:
        return IndexRange{1, nx, 1, ny};
    
    case Staggering::UFace:
        return IndexRange{1, nx-1, 1, ny};

    case Staggering::VFace:
        return IndexRange{1, nx, 1, ny-1};
    }

    assert(false && "unknown staggering");

    return IndexRange{};

}

void laplacian(const Field& f, const Grid& g, Staggering s, Field& out){

    assert((out.getNi() == f.getNi()) && (out.getNj() == f.getNj()));

    const auto [i0, i1, j0, j1] = interiorRange(g, s);
    
    const double dx = g.getDx();
    const double dy = g.getDy();

    // Precompute coefficients
    const double idx2 = 1.0 / (dx * dx);
    const double idy2 = 1.0 / (dy * dy);

    for (int j = j0; j <= j1; ++j){

        for (int i = i0; i <= i1; ++i){

            const double c = f(i, j);
            out(i, j) = (f(i+1, j) - 2.0*c + f(i-1, j)) * idx2
                        + (f(i, j+1) - 2.0*c + f(i, j-1)) * idy2;

        }

    }
}

double derivative(const Field& f, int i, int j, int di, int dj, double h, double a, AdvectionScheme scheme){

    // 1. Central difference: always available (1 neighbour each side)
    const double central = ( f(i+di, j+dj) - f(i-di, j-dj) ) / (2*h);

    if (scheme == AdvectionScheme::Central){
        return central;
    }

    // 2. Upwind2 needs 2 neighbours. Do they exist on BOTH sides?
    const bool fits = (i - 2*di >= 0) && (i + 2*di <= f.getNi() - 1)
                    && (j - 2*dj >= 0) && (j + 2*dj <= f.getNj() - 1);
    if (!fits){
        return central;                        // fallback near the boundary
    }
    // 3. Pick the upstream side by the sign of a
    if (a > 0){
        return (3 * f(i,j) - 4 * f(i-di, j-dj) + f(i-2*di, j-2*dj) ) / (2*h);
    }
    else{
        return ( -3 * f(i,j) + 4 * f(i+di, j+dj) - f(i+2*di, j+2*dj) ) / (2*h);
    }
}

void advectionU(const Field& u, const Field& v, const Grid& g, AdvectionScheme scheme, Field& out){

    assert((u.size() == v.size()) && (u.size() == out.size()));
    
    const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::UFace);
    for (int j = j0; j <= j1; ++j){
        for (int i = i0;  i<= i1; ++i){

            double a = u(i, j);
            double b = vAtU(v, i, j);
            out(i, j) = a * derivative(u, i,j, 1,0, g.getDx(), a, scheme) + b * derivative(u, i,j, 0,1, g.getDy(), b, scheme);

        }
    }

}
void advectionV(const Field& u, const Field& v, const Grid& g, AdvectionScheme scheme, Field& out) {

    assert((u.size() == v.size()) && (u.size() == out.size()));
    
    const auto [i0, i1, j0, j1] = interiorRange(g, Staggering::VFace);
    for (int j = j0; j <= j1; ++j){
        for (int i = i0; i <= i1; ++i){

            const double a = uAtV(u, i, j);
            const double b = v(i, j);
            out(i, j) = a * derivative(v, i,j, 1,0, g.getDx(), a, scheme) + b * derivative(v, i,j, 0,1, g.getDy(), b, scheme);

        }
    }
}


